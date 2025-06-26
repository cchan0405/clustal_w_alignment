# HIV-1 pol Gene Sequence Analysis

## Description  
This project analyzes a 1,302 bp partial coding sequence of the *HIV-1 pol* gene from a ViroSeq patient isolate (Cameroon1), focusing on regions associated with antiretroviral drug resistance. The analyzed fragment likely includes the complete **protease** domain and a portion of the **reverse transcriptase** domain.

## Features  
- DNA-to-protein translation using Biopython  
- Multiple sequence alignment with ClustalW  
- Alignment score extraction and visualization  
- Phylogenetic tree construction  
- Command-line interface via `argparse`  

## Usage  
All code is contained in `reading.py`. Input FASTA files (patient sequences and reference) should be in the same folder. Run the script as follows:

```bash
python reading.py --input sequence.fasta --ref reference_sequence.fasta
```
## To view available command-line options:

```bash
python reading.py --help
```
## Requirements
- Python 3.7 or higher
- Biopython
- Matplotlib
- ClustalW (Note: I installed ClustalW through Anaconda)

## Notes
- Output files include translated protein sequences (proteins.fasta), alignment (proteins.aln), phylogenetic tree file (proteins.dnd), and an alignment score plot (alignment_scores.png).
- ClustalW can be installed via conda:
```bash
conda install -c conda-forge -c bioconda clustalw
```
