from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio import AlignIO
import subprocess
import matplotlib.pyplot as plt
import argparse
import shutil
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Translate sequences, align with clustalw, and plot alignment scores. You can download clustalw off of anaconda, assuming you don't have it yet.  Make sure to run --input and --ref to using the parser")
    parser.add_argument('--input', default='/Users/chesterchan/Downloads/Python/HIV/sequence.fasta', help='Input FASTA file with DNA sequences')
    parser.add_argument('--ref', default='/Users/chesterchan/Downloads/Python/HIV/reference_sequence.fasta', help='Reference FASTA file')
    parser.add_argument('--out', default='proteins.fasta', help='Output FASTA file for translated proteins')
    args = parser.parse_args()

    if shutil.which("clustalw") is None:
        print("Error: ClustalW is not installed.")
        sys.exit(1)

    sequences = list(SeqIO.parse(args.input, "fasta"))

    reference_sequence = next(SeqIO.parse(args.ref, "fasta"))
    ref_seq_trimmed = reference_sequence.seq[:len(reference_sequence.seq) - len(reference_sequence.seq) % 3]
    reference_protein = ref_seq_trimmed.translate(to_stop=True)

    protein_records = [SeqRecord(reference_protein, id="reference", description="translated reference protein")]

    for record in sequences:
        seq_trimmed = record.seq[:len(record.seq) - len(record.seq) % 3]
        protein = seq_trimmed.translate(to_stop=True)
        protein_record = SeqRecord(protein, id=record.id, description="translated patient protein")
        protein_records.append(protein_record)

    SeqIO.write(protein_records, args.out, "fasta")

    cmd = ["clustalw", f"-infile={args.out}"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stderr:
        print("ClustalW STDERR:")
        print(result.stderr)

    scores = []
    for line in result.stdout.splitlines():
        if line.startswith("Sequences (1:") and "Aligned. Score:" in line:
            print(line.strip())
            score_str = line.strip().split("Score:")[-1].strip()
            try:
                scores.append(int(score_str))
            except ValueError:
                print(f"Could not convert score to int: {score_str}")
    tree_file = os.path.splitext(args.out)[0] + ".dnd"
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree, do_show=False)

    aln_file_read = os.path.splitext(args.out)[0] + ".aln"
    alignment = AlignIO.read(aln_file_read, "clustal")
    print(alignment)
    
    plt.figure(figsize=(12, 4))
    plt.plot(range(2, len(scores) + 2), scores, marker='o')
    plt.title("Alignment Scores: Reference vs. Other Sequences")
    plt.xlabel("Sequence Number")
    plt.ylabel("Alignment Score")
    plt.grid(True)

    plt.show()

    

if __name__ == "__main__":
    main()
