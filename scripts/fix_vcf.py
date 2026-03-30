#!/usr/bin/env python3
"""Resolve symbolic SV alleles in a VCF stream.

Replaces <DEL> with the reference sequence and drops unresolved symbolic
alleles (<INV>, <DUP>, <INS>). Reads from stdin, writes to stdout.

Usage: zcat input.vcf.gz | python fix_vcf.py <reference.fa> | bgzip > output.vcf.gz
"""
import sys
import pysam

ref_fn = sys.argv[1]
ref = pysam.FastaFile(ref_fn)

vcf = pysam.VariantFile("-")
out = pysam.VariantFile("-", 'w', header=vcf.header)

for entry in vcf:
    if not entry.alts:
        continue
    if entry.alts[0] == '<DEL>':
        entry.ref = ref.fetch(entry.chrom, entry.start, entry.stop)
        entry.alts = [entry.ref[0]]
    elif entry.alts[0].startswith('<'):
        continue
    out.write(entry)
