#!/bin/bash
echo "Setting up data"

# Sequencing cost dataset
SEQUENCING_DATA="https://www.genome.gov/sites/default/files/media/files/2021-11/Sequencing_Cost_Data_Table_Aug2021.xls"
if [ ! -f assets/downloads/Sequencing_Cost_Data_Table_Aug2021.xls ]; then
    mkdir -p assets/downloads
    wget -P assets/downloads ${SEQUENCING_DATA}
fi

TSKIT_DOWNLOAD="https://tskit.dev/tutorials/examples/download.py"
if [ ! -f data/basics.trees ]; then
    wget ${TSKIT_DOWNLOAD}
    python download.py
    rm -f download.py
fi
