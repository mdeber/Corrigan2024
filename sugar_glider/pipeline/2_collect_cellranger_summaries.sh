#!/usr/bin/env bash

cd $(dirname $0)/..
ALIGN_DIR=`pwd`/cellranger
[ -e metadata ] || mkdir metadata

# collect summary.csv files into a single csv file
touch metadata/cellranger_summaries.csv
cat samples.txt | while read SAMPLE; do
    PATH_CSV="$ALIGN_DIR/$SAMPLE/$SAMPLE/outs/metrics_summary.csv"
    [ -e $PATH_CSV ] || continue
    HEADER=$(head -1 $PATH_CSV)
    DATA=$(grep -v "$HEADER" $PATH_CSV)
    if [ ! -s metadata/cellranger_summaries.csv ]; then
        echo "Sample ID,$HEADER" >> metadata/cellranger_summaries.csv
    fi
    echo "$SAMPLE,$DATA" >> metadata/cellranger_summaries.csv
done
exit 0
