#!/usr/bin/env python
# 
# Make new tblastx mapping results that use gene names instead of transcript/CDS identifiers.
#
# For each gene pair in tblastx results, choose the transcript pair with the highest score;
# note that this script won't record which specific transcripts were chosen, it just returns
# the max blast score for each gene pair and writes out those new maps
# 
import pandas as pd
import numpy as np
from pathlib import Path
import os
import re

script_dir          = Path(__file__).parent.resolve()
dir_samap           = script_dir.parent.resolve()
pdir_maps           = Path(dir_samap, "maps")
pdir_genes_for_maps = Path(dir_samap, "genes_for_maps")
pdir_maps_out       = Path(dir_samap, "maps_by_gene")


# --------------------------------------------------
# Import the maps
# --------------------------------------------------

print("Importing original transcript maps...")
paths_maps = {
    re.sub(".*/(.*).txt$", "\\1", str(nm)) : nm for nm in pdir_maps.rglob("*/*")
}
# all_maps = {k: pd.read_csv(v, sep='\t', index_col=0, header=None) for k, v in paths_maps.items()}
all_maps = {k: pd.read_csv(v, sep='\t', header=None) for k, v in paths_maps.items()}

# >>> all_maps.keys()
# dict_keys(['hg_to_mm', 'mm_to_hg', 'hg_to_cj', 'cj_to_hg', 'cj_to_mm', 'mm_to_cj'])
# 
# >>> all_maps['hg_to_cj']
#                                                  0                                         1       2    3    4   5    6    7    8    9              10     11
# 0             lcl|NC_000001.11_cds_NP_001005484.2_1  lcl|NC_071451.1_cds_XP_009003880.3_44831  86.557  305   41   0   66  980    3  917   0.000000e+00  659.0
# 1             lcl|NC_000001.11_cds_NP_001005484.2_1  lcl|NC_071451.1_cds_XP_017832442.1_44839  59.450  291  118   0   76  948   34  906  1.900000e-117  422.0
# 2             lcl|NC_000001.11_cds_NP_001005484.2_1  lcl|NC_071451.1_cds_XP_017832440.1_44838  59.450  291  118   0   76  948   34  906  1.900000e-117  422.0
# 3             lcl|NC_000001.11_cds_NP_001005484.2_1  lcl|NC_071451.1_cds_XP_009004114.2_44837  59.450  291  118   0   76  948   34  906  1.900000e-117  422.0
# 4             lcl|NC_000001.11_cds_NP_001005484.2_1  lcl|NC_071451.1_cds_XP_017832439.1_44836  59.450  291  118   0   76  948   34  906  1.900000e-117  422.0
# ...                                             ...                                       ...     ...  ...  ...  ..  ...  ...  ...  ...            ...    ...
# 11696645     lcl|NT_113949.2_cds_NP_703144.3_123397  lcl|NC_071463.1_cds_XP_054106070.1_79608  36.364   66   42   0  283  480  289  486   2.720000e-16   51.9
# 11696646     lcl|NT_113949.2_cds_NP_703144.3_123397  lcl|NC_071463.1_cds_XP_054106077.1_79674  39.130   69   42   0  151  357  358  564   1.110000e-12   60.6
# 11696647     lcl|NT_113949.2_cds_NP_703144.3_123397  lcl|NC_071463.1_cds_XP_035141267.2_79675  40.678   59   35   0  619  795  319  495   2.930000e-13   55.1
# 11696648     lcl|NT_113949.2_cds_NP_703144.3_123397  lcl|NC_071463.1_cds_XP_002762527.1_79652  36.207   58   37   0   70  243   70  243   9.910000e-10   45.5
# 11696649  lcl|NC_012920.1_cds_YP_003024035.1_123407   lcl|NC_071442.1_cds_XP_054096692.1_5603  37.500   40   25   0  738  857  763  882   1.730000e-10   36.3

# [11696650 rows x 12 columns]

### this would be the "BLASTn tabular output format 6", wich I think might be [from tblastx -outfmt 6]
#   0.  qseqid      query or source (gene) sequence id
#   1.  sseqid      subject or target (reference genome) sequence id
#   2.  pident      percentage of identical positions
#   3.  length      alignment length (sequence overlap)
#   4.  mismatch    number of mismatches
#   5.  gapopen     number of gap openings
#   6.  qstart      start of alignment in query
#   7.  qend        end of alignment in query
#   8.  sstart      start of alignment in subject
#   9.  send        end of alignment in subject
#  10.  evalue      expect value
#  11.  bitscore    bit score

# --------------------------------------------------
# Import gene names for transcript elements as a mapping dictionary
# --------------------------------------------------

print("Importing the mapping between CDS/transcript names and gene names...")
genes_for_maps = {
    re.sub(".csv$", "", nm) : pd.read_csv(Path(pdir_genes_for_maps, nm), header=None, index_col=0).to_dict()[1] for nm in os.listdir(pdir_genes_for_maps)
}

# >>> genes_for_maps.keys()
# dict_keys(['mm', 'cj', 'hg'])

# >>> genes_for_maps['mm']['lcl|NC_000081.6_cds_NP_034577.1_70576']
# 'Hnrnpa1'

# --------------------------------------------------
# Map identifiers to gene names
# --------------------------------------------------

print("Mapping CDS/transcript names to their gene names...")

# I'm copying and leaving intermediates just to simplify finding the chosen transcripts later if needed
all_maps_genes = all_maps.copy()

for nm_pair in list(all_maps_genes.keys()):
    nm1 = re.sub("_.*", "", nm_pair)
    nm2 = re.sub(".*_", "", nm_pair)
    all_maps_genes[nm_pair][0] = all_maps_genes[nm_pair][0].map(genes_for_maps[nm1])
    all_maps_genes[nm_pair][1] = all_maps_genes[nm_pair][1].map(genes_for_maps[nm2])

# --------------------------------------------------
# For each gene pair, take mapping with max score
# --------------------------------------------------

print("For each gene pair, selecting the mapping with the highest tblastx score...")
final_maps_genes = all_maps_genes.copy()
for nm_pair in list(final_maps_genes.keys()):
    final_maps_genes[nm_pair] = final_maps_genes[nm_pair].sort_values(by=11, ascending=False).groupby([0,1], group_keys=False, sort=False).first()

# --------------------------------------------------
# Write out new maps, matching folder/file structure for SAMap
# --------------------------------------------------

print("Writing out the new gene-to-gene mappings...")
if not pdir_maps_out.is_dir():
    os.mkdir(pdir_maps_out)

for nm_pair in os.listdir(pdir_maps):
    if not Path(pdir_maps, nm_pair).is_dir():
        continue
    assert (len(nm_pair) == 4)
    if not Path(pdir_maps_out, nm_pair).is_dir():
        os.mkdir(Path(pdir_maps_out, nm_pair))
    nm1 = nm_pair[0:2]
    nm2 = nm_pair[2:4]
    for f in [f"{nm1}_to_{nm2}", f"{nm2}_to_{nm1}"]:
        out = Path(pdir_maps_out, nm_pair, f + ".txt")
        final_maps_genes[f].to_csv(out, sep="\t", header=False)

print("Script finished.")
