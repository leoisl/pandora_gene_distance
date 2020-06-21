import pandas as pd
df = pd.read_csv("/hps/nobackup/iqbal/leandro/pandora1_paper/metadata/variant_calls.20_way.pandora.snippy.samtools.illumina.100x.csv")
samtools_df = df[df.tool.apply(lambda tool: tool.startswith("samtools"))]
samtools_df = samtools_df[samtools_df.sample_id == "063_STEC"]
samtools_df["id"] = samtools_df["tool"].apply(lambda tool : tool.replace("samtools_", ""))
samtools_df["fasta"] = samtools_df["reference"]
samtools_df[["id", "fasta"]].to_csv("references_samtools.csv", index=False)