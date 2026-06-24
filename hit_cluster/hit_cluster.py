#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict
from Bio import SeqIO
import subprocess

##fazer argparse
def parse_args():
    parser = argparse.ArgumentParser(
        description="Join overlapping hits and avoid loci repetition."
    )
    parser.add_argument("Query", help="Path to the input fasta file (note that hits will be clustered independt of query sequence)")
    parser.add_argument("Subject", help="Path to DB files used in BLAST")
    parser.add_argument("-ex", type=int, default=0, help="expanding size (both sides), default = 0")
    parser.add_argument("-cy", type=int, default=0, help="expanding cycles, default = 0")
    parser.add_argument("-o", type=str, default="output.fa", help="output file (fasta)")
    parser.add_argument("-s", type=str, default="\t", help="type of separator used in the file, default = tab")

    return parser.parse_args()

##Run BLAST
def blast_run(db_file, query_file, step):
  tsv_file = "blast_results.tsv"
  subprocess.run(["makeblastdb", "-in", db_file, "-dbtype", "nucl"])

  if step == 1: ## Running blast to first ID loci
    try:
      subprocess.run(["blastn", "-query", query_file, "-db", db_file, "-outfmt", "6 qseqid sseqid sstart send evalue", "-out", tsv_file])
    except Exception as e:
      print(f"Error: {e}, first BLAST didn't run accordingly, examine files.")

  if step == 2: ##Running blast to fecth hit coordinates inside expanded and grouped loci
    try:
      step2_out = "hit_coords_in_loci.bed"
      subprocess.run(["blastn", "-query", query_file, "-db", db_file, "-outfmt", "6 sseqid sstart send qseqid", "-out", step2_out])
    except Exception as e:
      print(f"Error: {e}, final BLAST didn't run accordingly, examine files.")

##parsing tsv using file separator (e.g TSV -> \t)
def processing_tsv(tsv_file, sep):
    loci = set()
    hits = []
    groups = defaultdict(list)
    with open(tsv_file,"r") as f:
      for line in f:
          # Assuming format: Query-ID \sep Subject-id \sep Subject Start \sep Subject End \sep anything else
          parts = line.strip().split(f"{sep}")
          if len(parts) < 5:
              print(f"Skipping malformed line: {line.strip()} (expected at least 6 fields)")
              continue
          q = parts[0] # Query-ID
          # parts[1] is 'Subject-id'
          s = parts[1] # Subject scaffold number
          start = int(parts[2])
          end   = int(parts[3])
          e = float(parts[4]) # E-value

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
    if not groups[s]: # Skip if no hits for this subject
        continue

    # Initialize with the first hit for the current subject
    current_q, current_start, current_end, current_e = groups[s][0]

    for i in range(1, len(groups[s])):
        next_q, next_start, next_end, next_e = groups[s][i]
        # Check for overlap
        if overlaps(current_start, current_end, next_start, next_end) or overlaps(next_start, next_end, current_start, current_end):
            # Merge overlapping ranges: take min of starts, max of ends.
            # Keep q and e from the first hit that initiated the merged block.
            current_start = min(current_start, next_start)
            current_end = max(current_end, next_end)

        else:
            # If no overlap, add the current merged range to the clustered_groups
            clustered_groups[s].append((current_q, current_start, current_end, current_e))
            # Start a new merged range with the non-overlapping hit
            current_q, current_start, current_end, current_e = next_q, next_start, next_end, next_e

    # Add the last (or only) processed range for the current subject
    clustered_groups[s].append((current_q, current_start, current_end, current_e))

  return clustered_groups

def expand_ranges(groups, expanding_size):
  expanded_groups = defaultdict(list)
  for s in groups:
    for q, start, end, e in groups[s]:
      expanded_start = max(0, start - expanding_size)
      expanded_end = end + expanding_size
      # Preserve all original hit information along with expanded coordinates
      expanded_groups[s].append((q, expanded_start, expanded_end, e))
  return expanded_groups

def write_tsv(clustered_groups):
  with open("clustered_ranges.bed", "w") as f:
    # Iterate through each subject (s) and its list of clustered ranges
    for s, ranges_list in clustered_groups.items():
      # Iterate through each clustered range (q, start, end, e) for the current subject
      for q, start, end, e in ranges_list:
          f.write(f"{s}\t{start}\t{end}\t{q}\n")


def get_fasta(grouped_tsv_file, subject, out_file):
  subprocess.run(["bedtools", "getfasta", "-fi", subject, "-bed", grouped_tsv_file, "-fo", out_file])
  return out_file

def main():
  #parsing argv
  args = parse_args()
  query_file = args.Query
  db_file = args.Subject
  blast_run(db_file, query_file, 1)
  tsv_file = "blast_results.tsv"
  expanding_cycles = args.cy
  expanding_size = args.ex
  out_file = args.o
  sep = args.s

  # Initial processing of the TSV file
  groups = processing_tsv(tsv_file, sep)
  groups = compare_ranges(groups)

  # Perform iterative expansion and clustering if specified

  while expanding_cycles > 0:
    expanding_cycles -= 1
    expanded_groups = expand_ranges(groups, expanding_size)
    groups = compare_ranges(expanded_groups) # Update 'groups' with the clustered results for the next iteration

  # Write the final clustered ranges to a BED file

  write_tsv(groups)
  get_fasta("clustered_ranges.bed", db_file, out_file)
  blast_run(out_file, query_file, 2)


if __name__ == "__main__":
  main()
