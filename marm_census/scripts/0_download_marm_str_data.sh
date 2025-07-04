MARM_SCRIPT_DIR=$(dirname $0);
MARM_DIR=${MARM_SCRIPT_DIR}/..;

[ -e ${MARM_DIR}/data ] || mkdir ${MARM_DIR}/data;

cd ${MARM_DIR}/data;

# 
# Data from manuscript:
#   Fenna M. Krienen et al. 
#   A marmoset brain cell census reveals regional specialization of cellular identities.
#   Sci. Adv.9,eadk3986(2023). DOI:10.1126/sciadv.adk3986
# 
# But as of writing, the only publically available processed data (posted to CZI CellxGene)
# contains only genes that have 1:1 orthologs in humans, and so is not the original dataset
# used in the publication:
#   https://cellxgene.cziscience.com/collections/0fd39ad7-5d2d-41c2-bda0-c55bde614bdb
#
# Note CZI accepted a non-humanized version of the data in 2025, but that version uses
# ENSEMBL instead of RefSeq, which is far under-annotated for marmoset, so is still not
# the same as what we analyzed.
#
# While we find a solution to this, the original "striatum.GAD" dataset from that publication
# (converted to Seurat RDS) will be included in this GitHub repo.
#
exit 0
