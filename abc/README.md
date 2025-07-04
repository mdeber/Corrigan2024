The entire ABC atlas can be cloned from `arn:aws:s3:::allen-brain-cell-atlas`
(i.e. https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html)

Most analyses only need these files:
`allen-brain-cell-atlas/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-STR-raw.h5ad`
https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-STR-raw.h5ad 

`allen-brain-cell-atlas/metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv`
https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv

But a number of other datasets are also necessary for the final scVI integration, which included additional cells for comparative purposes. Run script `abc/scripts/0_download_abc_str_data.sh` to download all files into `abc/data`.

(For initial analyses, it would have saved time & space to subset to the class 08 cells in python before making the seurat rds object).
