Corrigan, DeBerardine et al. 2025

The main adult comparative analyses are in scripts in the top-level directory, except for `SAMap`, which contains its own pipeline with analysis.
Additionally, the preprocessing of the `sugar_glider` data is here, but isn't a part of the comparative analyses.

Before running `scvi_integrations.py`, `barplot_celltype_prevalence.R`, or `correlation_analysis.R`, all scripts in `abc`, `marm_census`, and `siletti` must be run.
See details about `marm_hmba` in that folder. Cell barcodes and annotations are provided here for that dataset. The presence of those annotations & cell barcodes means everything besides the scVI integration can be run without downloading/processing that data.

Running `reciprocal_hotspot_projections.R` requires that `abc`, `marm_census`, `siletti`, `dev_ferret`, `dev_macaque`, `dev_mouse`, and `dev_pig` all have their pipelines' run.

If run from the command line (e.g., using `bash <script>.sh`, `Rscript <script>.R`, and `python <script.py>`), all scripts will automatically find the local paths of all required files, provided you don't reorganize any files in this repo. But beware, most of the bash scripts might not find their locations if launched in slurm jobs (you should substitute in `realpath` in those cases case).

The preprint used liger integrations, which we eventually decided against for the sake of brevity and because scVI performed better, despite our wide exploration of liger parameters.
