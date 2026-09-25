"""Clean Alena's lifted-over region set (hg -> mm10, UCSC web liftOver output) before counting.

Steps:
  1. drop multi-mapped inputs entirely. With "allow multiple output regions", UCSC writes
     every hit of one input as consecutive rows scored 1, 2, 3...; the hits are not
     ranked, so hit 1 is not necessarily the right one (--keep-multi keeps all hits)
  2. drop chrM and exact duplicates
  3. optional width filter on the lifted widths (--min-width / --max-width), for use
     once Alena confirms the original human widths and liftOver minMatch
  4. resize every region to a fixed width around its midpoint
  5. merge regions that overlap after resizing, so no fragment is counted twice
  6. optionally drop regions overlapping the ENCODE mm10 blacklist (--blacklist)

Usage:
  python prep_regions.py Table_S13_hglft_genome_248500_fceb0.bed regions_mm10_clean.bed \
      --width 500 --blacklist mm10-blacklist.v2.bed
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("bed_in")
ap.add_argument("bed_out")
ap.add_argument("--width", type=int, default=500)
ap.add_argument("--keep-multi", action="store_true", help="keep all hits of multi-mapped inputs")
ap.add_argument("--min-width", type=int, default=None, help="drop lifted regions narrower than this")
ap.add_argument("--max-width", type=int, default=None, help="drop lifted regions wider than this")
ap.add_argument("--blacklist", help="BED of blacklisted regions to remove (e.g. ENCODE mm10 v2)")
args = ap.parse_args()

cols = ["chr", "start", "end", "name", "score", "strand"]
df = pd.read_csv(args.bed_in, sep="\t", header=None, names=cols, comment="#")
w0 = df.end - df.start
print(f"input: {len(df)} rows; width median {w0.median():.0f} bp, "
      f"<50 bp: {(w0 < 50).sum()}, >2 kb: {(w0 > 2000).sum()}, max {w0.max()}")

# each input region starts a new group at score 1; later hits (2, 3...) follow directly
group = (df.score == 1).cumsum()
gsize = group.map(group.value_counts())
n_multi_inputs = (gsize[df.score == 1] > 1).sum()
if not args.keep_multi:
    df = df[gsize == 1]
print(f"multi-mapped inputs: {n_multi_inputs} ({'kept' if args.keep_multi else 'dropped, all hits'}); "
      f"{len(df)} rows left")

df = df[df.chr != "chrM"].drop_duplicates(["chr", "start", "end"])
w = df.end - df.start
if args.min_width is not None:
    df = df[w >= args.min_width]
if args.max_width is not None:
    df = df[(df.end - df.start) <= args.max_width]
print(f"after chrM/duplicate/width filters: {len(df)}")

mid = (df.start + df.end) // 2
df = df.assign(start=(mid - args.width // 2).clip(lower=0))
df = df.assign(end=df.start + args.width).sort_values(["chr", "start"])

merged = []
for chrom, g in df.groupby("chr", sort=True):
    cs = ce = None
    for s, e in zip(g.start, g.end):
        if cs is None:
            cs, ce = s, e
        elif s <= ce:
            ce = max(ce, e)
        else:
            merged.append((chrom, cs, ce))
            cs, ce = s, e
    merged.append((chrom, cs, ce))
out = pd.DataFrame(merged, columns=["chr", "start", "end"])
print(f"after resizing to {args.width} bp and merging overlaps: {len(out)}")

if args.blacklist:
    bl = pd.read_csv(args.blacklist, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chr", "start", "end"])
    hit = pd.Series(False, index=out.index)
    for chrom, g in out.groupby("chr"):
        b = bl[bl.chr == chrom]
        for i, s, e in zip(g.index, g.start, g.end):
            hit.loc[i] = bool(((b.start < e) & (b.end > s)).any())   # half-open BED overlap
    out = out[~hit].reset_index(drop=True)
    print(f"removed {int(hit.sum())} regions overlapping the blacklist")

out["name"] = [f"reg{i:05d}" for i in range(1, len(out) + 1)]
out.to_csv(args.bed_out, sep="\t", header=False, index=False)
print(f"output: {len(out)} regions ({(out.end - out.start > args.width).sum()} longer than "
      f"{args.width} bp after merging) -> {args.bed_out}")
print("chromosomes:", out.chr.value_counts().to_dict())
