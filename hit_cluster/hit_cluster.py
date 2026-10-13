#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Join overlapping hits and avoid loci repetition while respecting FASTA boundaries."
    )
    parser.add_argument(
        "tsv_file",
        help="Path to the input tsv file (fmt = Query-ID \t Subject-id \t Subject \t Start \t Subject End \t anything else)",
    )
    parser.add_argument(
        "-f", "--fasta", required=True, help="Path to the reference FASTA file to check max sequence lengths."
    )
    parser.add_argument(
        "-ex", type=int, default=0, help="expanding size (both sides)"
    )
    parser.add_argument("-cy", type=int, default=0, help="expanding cycles")

    return parser.parse_args()


def get_sequence_lengths(fasta_file):
    """Reads the FASTA file and returns a dictionary of {sequence_id: length}."""
    seq_lengths = {}
    print(f"Reading sequence lengths from {fasta_file}...")
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths[record.id] = len(record)
    return seq_lengths


def processing_tsv(tsv_file):
    loci = set()
    hits = []
    groups = defaultdict(list)
    with open(tsv_file, "r") as f:
        for line in f:
            # Assuming format: Query-ID \t Subject-id \t Subject \t Start \t Subject End \t anything else
            parts = line.strip().split("\t")
            if len(parts) < 6:
                print(
                    f"Skipping malformed line: {line.strip()} (expected at least 6 fields)"
                )
                continue
            q = parts[0]  # Query-ID
            s = parts[2]  # Subject (locus identifier)
            start = int(parts[3])
            end = int(parts[4])
            e = float(parts[5])  # e-value

            # normalize orientation
            if start > end:
                start, end = end, start

            hits.append((q, s, start, end, e))
            loci.add(s)

        # Re-grouping all hits by subject after initial parsing and normalization
        for q, s, start, end, e in hits:
            groups[s].append((q, start, end, e))

        # Sort hits within each subject by start coordinate
        for s in groups:
            groups[s].sort(key=lambda x: x[1])
    return groups


def overlaps(a_start, a_end, b_start, b_end):
    if b_start <= a_end and b_end >= a_start:
        return True
    else:
        return False


def compare_ranges(groups):
    clustered_groups = defaultdict(list)
    for s in groups:
        if not groups[s]:  
            continue

        current_q, current_start, current_end, current_e = groups[s][0]

        for i in range(1, len(groups[s])):
            next_q, next_start, next_end, next_e = groups[s][i]

            if overlaps(current_start, current_end, next_start, next_end):
                current_start = min(current_start, next_start)
                current_end = max(current_end, next_end)
            else:
                clustered_groups[s].append(
                    (current_q, current_start, current_end, current_e)
                )
                current_q, current_start, current_end, current_e = (
                    next_q,
                    next_start,
                    next_end,
                    next_e,
                )

        clustered_groups[s].append(
            (current_q, current_start, current_end, current_e)
        )

    return clustered_groups


def expand_ranges(groups, expanding_size, seq_lengths):
    expanded_groups = defaultdict(list)
    for s in groups:
        # Get the max length for this chromosome/subject, default to infinity if not found in FASTA
        max_len = seq_lengths.get(s, float('inf'))
        
        for q, start, end, e in groups[s]:
            expanded_start = start - expanding_size
            if expanded_start < 0:
                expanded_start = 0
                
            expanded_end = end + expanding_size
            # Prevent the end coordinate from exceeding the actual sequence length
            if expanded_end > max_len:
                expanded_end = max_len
                
            expanded_groups[s].append((q, expanded_start, expanded_end, e))
    return expanded_groups


def write_tsv(clustered_groups):
    with open("clustered_ranges.tsv", "w") as f:
        for s, ranges_list in clustered_groups.items():
            for q, start, end, e in ranges_list:
                f.write(f"{s}\t{q}\t{start}\t{end}\n")


def main():
    args = parse_args()
    tsv_file = args.tsv_file
    fasta_file = args.fasta
    expanding_cycles = args.cy
    expanding_size = args.ex

    # Get maximum coordinates from FASTA
    seq_lengths = get_sequence_lengths(fasta_file)

    groups = processing_tsv(tsv_file)

    while expanding_cycles > 0:
        expanding_cycles -= 1
        # Pass the seq_lengths dictionary to the expand_ranges function
        expanded_groups = expand_ranges(groups, expanding_size, seq_lengths)
        groups = compare_ranges(expanded_groups)

    write_tsv(groups)


if __name__ == "__main__":
    main()
