#!/usr/bin/env python3
"""Convert 2bp CpG bedGraph coordinates to single-base cytosine coordinates.

Requires a FASTA file (from bedtools getfasta) to determine strand orientation.

Usage: python bedconvert.py <bedgraph> <fasta> <output>
Example: python bedconvert.py 11_modkit/EMSeq_HG002_LAB01_REP01.bedGraph \
             11_modkit/EMSeq_HG002_LAB01_REP01.fa \
             11_modkit/EMSeq_HG002_LAB01_REP01.converted.bedGraph
"""
import sys

bedfile = sys.argv[1] if len(sys.argv) > 1 else sys.exit(f"Usage: {sys.argv[0]} <bedgraph> <fasta> <output>")
fastafile = sys.argv[2]
outfile = sys.argv[3]

# parse FASTA: map "chr:start-end" -> sequence
fasta_dict = {}
header = None
with open(fastafile) as ff:
    for line in ff:
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            if "::" in header:
                header = header.split("::", 1)[-1]
            fasta_dict[header] = []
        else:
            fasta_dict[header].append(line.upper())

# parse bedGraph and convert coordinates
with open(bedfile) as bf, open(outfile, "w") as out:
    for line in bf:
        line = line.strip()
        if not line:
            continue
        cols = line.split("\t")
        chrom, start, end = cols[0], int(cols[1]), int(cols[2])
        key = f"{chrom}:{start}-{end}"

        if key not in fasta_dict:
            continue

        seq = "".join(fasta_dict[key])
        if seq == "CG":
            new_start, new_end = start, start + 1
        elif seq == "GC":
            new_start, new_end = start + 1, start + 2
        else:
            continue

        out.write("\t".join([chrom, str(new_start), str(new_end)] + cols[3:]) + "\n")

print(f"Converted: {outfile}")
