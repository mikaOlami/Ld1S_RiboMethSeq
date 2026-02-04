#!/usr/bin/env python
import os
import pandas as pd
from Bio import SeqIO
import numpy as np

# -------------------------------
# Resolve repo root from script location
# -------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Set working directory (where *.init / *.3p live)
directory = "."   # current working directory
ResultsDir = directory

# -------------------------------
# rRNA Nm annotation file (repo DB)
# -------------------------------
annot_file = os.path.join(
    REPO_DIR,
    "DB",
    "rRNA",
    "Ld1S_rRNA_annot_Nm.txt"
)
annot = pd.read_csv(annot_file, sep="\t")

# -------------------------------
# rRNA reference fasta file (repo DB)
# -------------------------------
fasta_path = os.path.join(
    REPO_DIR,
    "DB",
    "rRNA",
    "Ld1S_rRNA.fa"
)
myfasta = SeqIO.parse(fasta_path, "fasta")
rRNA_sequence = next(myfasta).seq

# -------------------------------
# Parameters
# -------------------------------
win_size = 6
rRNA_length = 8100

# -------------------------------
# Process each library
# -------------------------------
mylibs = [file.replace(".init", "") for file in os.listdir(directory) if file.endswith(".init")]

for thisLib in mylibs:
    # Read the table as it came from the pre-processing, must include the counts at the third column
    pre_myinit = np.loadtxt(thisLib + ".init", usecols=[2])
    pre_myends = np.loadtxt(thisLib + ".3p", usecols=[2])

    # Shifting reads 1 bp
    myinit = np.concatenate(([0], pre_myinit, np.zeros(rRNA_length - len(pre_myinit) - 1)))
    myends = np.concatenate((pre_myends[1:], np.zeros(rRNA_length - len(pre_myends) + 1)))
    mycov = myinit + myends
    mylength = len(myinit)

    Sc = np.full(mylength, np.nan)
    Sa = np.full(mylength, np.nan)
    Sb = np.full(mylength, np.nan)

    # Sum the weights for each position - avoid overweighted results
    W = (1 + (1 - 0.1 * win_size)) * win_size / 2

    # Calculate the scores
    for i in range(win_size, mylength - win_size - 1):
        # A score
        M_l = np.mean(mycov[i - win_size // 2:i])
        S_l = np.std(mycov[i - win_size // 2:i], ddof=1)
        M_r = np.mean(mycov[i + 1:i + win_size // 2 + 1])
        S_r = np.std(mycov[i + 1:i + win_size // 2 + 1], ddof=1)

        Sa[i] = max(
            0,
            1 - (2 * mycov[i] + 1)
            / (0.5 * abs(M_l - S_l) + mycov[i] + 0.5 * abs(M_r - S_r) + 1)
        )

        # B + C scores
        S1 = sum((1 - 0.1 * (j - 1)) * mycov[i - j] for j in range(1, win_size + 1))
        S1 /= W

        S2 = sum((1 - 0.1 * (j - 1)) * mycov[i + j] for j in range(1, win_size + 1))
        S2 /= W

        if S1 == 0 or S2 == 0:
            Sc[i] = 0
        else:
            Sc[i] = max(0, 1 - 2 * mycov[i] / (S1 + S2))

        Sb[i] = abs((mycov[i] - 0.5 * (S1 + S2)) / (mycov[i] + 1))

    # Extracting fasta base by base to be printed within the final table
    bp = list(str(rRNA_sequence)) + [np.nan] * (rRNA_length - len(rRNA_sequence))

    mydf = pd.DataFrame({
        "5p": myinit,
        "3p": myends,
        "cov": mycov,
        "Sa": Sa,
        "Sb": Sb,
        "Sc": Sc,
        "bp": bp
    })

    # Add annotation information
    mydf = pd.concat([mydf, annot], axis=1)

    # Write the final CSV file
    mydf.to_csv(os.path.join(ResultsDir, thisLib + ".csv"), index=False)
