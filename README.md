# Corrigan2024

The main analyses are done in `liger_integrations.R` and `reciprocal_hotspot_projections.R`. 

Before running `liger_integrations.R`, all scripts in `abc`, `marm_census`, and `siletti` must be run.
Before running `reciprocal_hotspot_projections.R`, all of those scripts must be run, as well as all scripts in `dev_ferret`, `dev_macaque`, `dev_mouse`, and `dev_pig`.

If run from the command line (e.g., using `bash <script>.sh`, `Rscript <script>.R`, and `python <script.py>`), all scripts will automatically find the local paths of all required files, provided you don't reorganize any files in this repo. 

Be aware the `liger_integrations.R` script is particularly poorly optimized and uses dozens of Gb of RAM and saves large files.
