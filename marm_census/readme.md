Note that the marmoset census striatum GABA neuron data is being included in this repo (see `scripts/0_download_marm_str_data.sh`).

The file `data/striatum.GAD.rds` was saved with LZMA (`saveRDS` with `compression="xz"`) which allows it to be hosted on GitHub but makes loading the file significantly slower. To make import faster, you can resave that file using `compression=FALSE`.
