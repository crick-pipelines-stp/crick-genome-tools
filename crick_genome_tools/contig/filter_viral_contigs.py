import pandas as pd
from itertools import combinations



def filter_viral_contigs(min_depth, min_ratio, contig_grouping, avdepth_path):
    contig_grouping = contig_grouping.split(',')
    df = pd.read_csv(avdepth_path, sep="\t")

    # Strip grouping from contig names to get base_id
    def strip_group(contig):
        s = contig
        for g in contig_grouping:
            s = s.replace(g, "")
        return s
    df["base_id"] = df["contig"].apply(strip_group)

    # Move base_id after contig
    cols = df.columns.tolist()
    cols.insert(1, cols.pop(cols.index("base_id")))
    df = df[cols]

    # --- Compute segment ratios within each base_id ---
    def find_group(contig):
        for g in contig_grouping:
            if g in contig:
                return g
        return "other"
    df["segment"] = df["contig"].apply(find_group)
    df["seg_ratio"] = df.groupby("segment")["av_depth"].transform(lambda x: x / x.max())
    df["seg_ratio"] = df["seg_ratio"].round(2)

    # --- Compute total depth and ratio to dominant base_id ---
    base_sums = df.groupby("base_id")["av_depth"].sum()
    max_depth = base_sums.max()
    base_ratios = (base_sums / max_depth).round(2)
    base_ratios.name = "ratio_to_max"
    df = df.merge(base_ratios, on="base_id")

    # Write data
    df.to_csv("contig_data.tsv", sep="\t", index=False)

    # --- Filter by min_depth ---
    df = df[df["av_depth"] >= min_depth]
    df = df[df["seg_ratio"] >= min_ratio]
    df = df[df["ratio_to_max"] >= min_ratio]

    valid_refs = df["contig"].tolist()
    with open("valid_refs.txt", "w") as f:
        f.write("\n".join(valid_refs))

