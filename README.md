# Camera: Constraint-based Assignment of Methyl Resonances

### Software dependencies

- `Python 3.4+` (additional modules below)
  - `NumPy`
  - `networkx`
  - `Pandas`
  - `halo`
  - `tqdm`
  - `biopython`

- `cryptominisat5` (`bin/env` must be able to find this binary)

### Preparing inputs

This program is used for assigning methyl resonances measured via NMR.
The files required to use this suite of programs are:
1. A crystal structure in PDB format (plus zero or more relaxations)
2. A 2D HMQC peak list in Sparky format
3. A 3D CCH peak list in Sparky format

The first step is to take these files and convert them into a format that is consumable by resonance assignment code.

1. PDB files:
the script `camera/predict_noes` will take as input one or more PDB files and output a `.json` file.
The file that the script outputs will contain a list of all methyls encountered in the PDB files the pairwise distances between all methyls in every structure.

2. HMQC peaklist:
the script `camera/scripts/hmqc2csv.py` takes as input a Sparky HMQC peaklist and outputs a `.csv` file containing the same information.
This suite uses `.csv` files to store NMR data for ease of parsing and easier human interactivity.

3. 3D CCH peaklist:
the script `camera/scripts/sparky2csv.py` takes as input a Sparky 3D CCH peaklist and outputs a `.csv` file containing the same information.
Make sure that the sparky list you export as input to the script includes 5 columns: a label, the three chemical shifts, and the intensity.

After performing these steps, you should have a `.json` file containing all information about the known structure, a `.csv` file describing all 2D peaks, and a `.csv` file describing all 3D peaks.

In some cases, you will have a 300ms and 50ms CCH experiment.
In this case, you want to mark the NOEs which occur in the 50ms experiment as "short" in the 3D `.csv` file.
To do this, run `camera/scripts/sparky2csv.py` on the 300ms and 50ms peak lists separately to create two `.csv` files, called `300ms.csv` and `50ms.csv`.
Executing the command `camera/mark_short_noes 300ms.csv 50ms.csv -o merged.csv` will merge the two files and indicate which NOEs were both observed and symmetrized in the 50ms experiment using a new column entitled `short`.

### Running symmetrization

The first step in assigning the methyl resonances is to process the NOEs.
The processing amounts to determining which NOEs are symmetric with each other in order to generate constraints.
The script will take as input the structure `.json` file emitted by `camera/predict_noes`, the 2D `.csv`, and the 3D `.csv` file.
The output will be a new 3D `.csv` file with a column labeled `reciprocal` which gives all possible reciprocals for each NOE.
Issuing the command `camera/process_noes hmqc.csv cch.csv structure.json -o processed.csv` will emit a list of processed NOEs with reciprocals labeled to `processed.csv`.

### Assigning the HMQC peaks

To run the assignment protocol, use the script `assign_signatures`.
This will take as input the structure `.json` file emitted by `camera/predict_noes`, the 2D `.csv` file, and the processed 3D `.csv` emitted by `camera/process_noes`.
The script will emit a 2D list with all possible assignments added for each 2D peak.
Issuing the command `camera/assign_signatures hmqc.csv processed.csv structure.json -o assigned_hmqc.csv` will generate a new 2D `.csv` file called `assigned_hmqc.csv` with an assignments column with the assignments for each peak.
