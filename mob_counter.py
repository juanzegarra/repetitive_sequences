#!/usr/bin/env python3
import sys
import os
import subprocess
from collections import defaultdict
from Bio import SeqIO

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file.fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    def mob_counter(fasta_file):
        mob_counter = defaultdict(int)
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.id)
            family = sequence.split("#")
            correct_family = family[1]

            if correct_family.startswith("DNA/"):
                if correct_family.startswith("DNA/hAT"):
                    correct_family = "hAT"
                else:
                    correct_family = correct_family.replace("DNA/", "")
                mob_counter[correct_family] += 1

            elif correct_family.startswith("PLE/"):
                correct_family = "PLE"
                mob_counter[correct_family] += 1

            elif correct_family.startswith("LTR/"):
                correct_family = "LTR"
                mob_counter[correct_family] += 1

            elif correct_family.startswith("LINE/"):
                correct_family = "LINE"
                mob_counter[correct_family] += 1

            elif correct_family.startswith("SINE/"):
                correct_family = "SINE"
                mob_counter[correct_family] += 1

            elif correct_family.startswith("RC/Helitron"):
                correct_family = "Helitron"
                mob_counter[correct_family] += 1

            elif correct_family.startswith("Satellite"):
                continue

            elif correct_family.startswith("simple_repeat"):
                continue

        file_name = os.path.basename(fasta_file).replace(".fasta", "")
        input_name = f"transposons_{file_name}.tsv"
        out_name = f"transposons_{file_name}.png"

        with open(input_name, "w") as file:
            file.write("Family\tCopies\n")
            for seq, count in mob_counter.items():
                file.write(f"{seq}\t{count}\n")

        subprocess.run(
            ["mob_counter_plot.R", input_name, out_name],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

    mob_counter(fasta_file)
    print(f"Processing {fasta_file}... done!")

if __name__ == "__main__":
    main()
