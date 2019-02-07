#!/bin/bash

# Print out all commands
set -x

# Convert 2D list from sparky to csv format
../scripts/hmqc2csv.py hmqc.list 2d.csv

# Convert 300ms CCH from sparky to CSV
../scripts/sparky2csv.py first_attempt_300ms.list 300ms.csv

# Convert 50ms CCH from sparky to CSV
../scripts/sparky2csv.py first_attempt_50ms.list 50ms.csv

# Extract a .json file form of the pdb entires
../predict_noes a01/*.pdb --chain A --colors ALA ILE LEU VAL -o g10.json

# Merge the NOE files
../mark_short_noes 300ms.csv 50ms.csv -o merged.csv

# Run symetrization
../process_noes 2d.csv merged.csv g10.json -o processed.csv

# Run assignment
../assign_signatures 2d.csv processed.csv g10.json -o assigned_2d.csv

