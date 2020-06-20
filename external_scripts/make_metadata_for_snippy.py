base_path="/hps/nobackup/iqbal/leandro/snippy_analysis_pipeline/processed_data/illumina/100x/random/063_STEC"

from glob import glob
import pandas as pd
from pathlib import Path

ids=[]
fastas=[]

for file in glob(f"{base_path}/*.ref.fa"):
    id = Path(file).name
    id = id.replace("snippy_063_STEC_AND_", "").replace(".ref.fa", "")
    ids.append(id)
    fastas.append(str(file))

df = pd.DataFrame({"id": ids, "fasta": fastas})
df.to_csv("references_snippy.csv", index=False)