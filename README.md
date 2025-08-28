# HIV Sequence Analysis and ClustalW Alignment

## Overview
This project provides Python scripts to **read, process, and align HIV sequences** using Biopython. It allows sequence parsing, multiple sequence alignment via **ClustalW**, and basic phylogenetic analysis.

---

## Features
- Parse DNA/RNA/Protein sequences from FASTA files using Biopython (`SeqIO`, `SeqRecord`).
- Perform multiple sequence alignment with **ClustalW**.
- Generate phylogenetic trees from aligned sequences using Biopython’s `Phylo` module.
- Save alignment outputs (`.aln` and `.dnd`) for further analysis.

---

## Requirements
- Python 3.11+
- Biopython
- ClustalW executable installed and accessible in your system PATH (Conda was used for this project)
