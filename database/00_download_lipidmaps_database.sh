#!/bin/bash

wget -O "raw_data/lipidmaps" "http://www.lipidmaps.org/data/structure/LMSDSearch.php?Mode=ProcessStrSearch&OutputMode=File&OutputQuote=No&OutputType=TSV&OutputColumnHeader=Yes"
python 00_preprocess_lipidmaps_database.py raw_data/lipidmaps raw_data/lipidmaps.tsv