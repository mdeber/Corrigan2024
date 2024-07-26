import pandas as pd
from pathlib import Path

path_ortho_tables = Path(__file__).parent.resolve()
path_df_ortho = Path(path_ortho_tables, "marm_human_mouse_ortholog_table_ens110.csv")
path_marm_id_table = Path(path_ortho_tables, "marm_ens110_ncbi_gene_id_conversion.txt")
path_out = Path(path_ortho_tables, "marm_mouse_orthologs.csv")

# ---
# Import gene name and ortholog lookup tables
# ---

marm_id_table = pd.read_csv(path_marm_id_table) # n=41577
marm_id_table.set_index('Gene stable ID', inplace = True)

df_ortho = pd.read_csv(path_df_ortho) # n=47028
df_ortho.set_index('Gene stable ID', inplace = True)

# ---
# Add marmoset gene names to the ortholog table
# ---

marm_id_table = marm_id_table.loc[marm_id_table['Gene name'].notna()] # n=25685
# from marm_id_table, only need the Gene name
marm_id_table = marm_id_table.loc[:, "Gene name"]

df_ortho = df_ortho.join(marm_id_table, how = "inner")
# sum(df_ortho.index.duplicated()) # n=2192

# ---
# Select orthologs
# ---

# Keep only high quality (n=16802 for humans; but only n=13084 for mouse)
df_ortho = df_ortho[
    (df_ortho['Mouse orthology confidence [0 low, 1 high]'] == 1)
]

# Then, to choose which ortholog to use:
# Easiest is to sort the dataframe by multiple columns according to our top priorities for the ortholog,
# then index by gene name, then keep only the "first" of each duplicated index (the highest priority ortholog)
df_ortho = df_ortho.sort_values(
    by = ['Gene name', 'Gene stable ID version', '%id. query gene identical to target Mouse gene'],
    ascending = [True, False, False]
)
df_ortho = df_ortho.set_index('Gene name')
df_ortho = df_ortho.loc[[not i for i in df_ortho.index.duplicated(keep='first')], ] # n=12671

# ---
# Write out
# ---

pd.DataFrame.to_csv(df_ortho, path_out)
