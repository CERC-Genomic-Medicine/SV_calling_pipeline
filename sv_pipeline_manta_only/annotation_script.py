#!/usr/bin/env python3
# last updated nov 2024 
import argparse
import numpy as np
from intervaltree import Interval, IntervalTree
import interval
import pysam
import os
import pandas as pd


def find_overlaps(bqc_filename, gnomad_filename, svtype):
    print(f"Processing BQC file: {bqc_filename}")
    print(f"Processing GnomAD file: {gnomad_filename}")
    print(f"SV type: {svtype}")

    overlaps = []  
    with pysam.VariantFile(bqc_filename) as bqc_vcf, pysam.VariantFile(gnomad_filename) as gnomad_vcf: 
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
        tree_bqc = IntervalTree()
        tree_gnomad = IntervalTree()
        for chrom in chromosomes:#by chromosome
            print(f"Processing: {chrom}")
            # Load all GnomAD SVs into interval tree
            for gnomad_variant in gnomad_vcf.fetch(chrom):
                gnomad_var_info = gnomad_variant.info
                if gnomad_var_info['SVTYPE'] == svtype:
                    if svtype == 'INS':
                        tree_gnomad.addi(gnomad_variant.pos, gnomad_variant.pos + gnomad_var_info['SVLEN'], data = gnomad_variant.info)
                    elif svtype == 'BND':
                        tree_gnomad.addi(gnomad_variant.pos, gnomad_variant.pos + 1, data = gnomad_variant.info)
                    else:
                        tree_gnomad.addi(gnomad_variant.pos, gnomad_variant.stop, data = gnomad_variant.info)
            # Iterate over study SVs
            bqc_iter = bqc_vcf.fetch(chrom)
            for bqc_variant in bqc_iter: #add each variant in the file to the interval tree
                # if bqc_variant.info['SVTYPE'] != svtype:
                #     continue
                start = bqc_variant.pos
                if svtype == 'INS': 
                    stop = bqc_variant.pos + bqc_variant.info['SVLEN']
                elif svtype == 'BND':
                    stop = bqc_variant.pos + 1 
                else:    
                    stop = bqc_variant.stop
                for interval in tree_gnomad.overlap(start, stop):
                    overlap_start = max(interval.begin, start)
                    overlap_end = min(interval.end, stop)
                    overlap_len = overlap_end - overlap_start
                    assert overlap_len > 0
                    if svtype == 'BND':
                        if (bqc_variant.info['CHR2'] == interval.data['CHR2']) and (bqc_variant.stop == interval.data['END2']):
                            ratio1 = 1.0
                            ratio2 = 1.0
                        else:
                            ratio1 = 0
                            ratio2 = 0
                    else:
                        ratio1 = overlap_len / abs(bqc_variant.info['SVLEN'])
                        ratio2 = overlap_len / interval.data['SVLEN']
                    if ratio1 >= 0.5 and ratio2 >= 0.5:
                        overlaps.append([svtype, 
                                         bqc_variant.id, 
                                         bqc_variant.chrom, 
                                         bqc_variant.pos, 
                                         bqc_variant.stop, 
                                         bqc_variant.info['AC'][0], 
                                         bqc_variant.info['AN'],
                                         bqc_variant.info['AC'][0]/bqc_variant.info['AN'],
                                         bqc_variant.info['SVLEN'],
                                         bqc_variant.info['F_MISSING'][0],
                                         interval.data['AC'][0], 
                                         interval.data['AN'],
                                         interval.data['AC'][0]/interval.data['AN'],
                                         interval.data['SVLEN'],
                                         ratio1,
                                         ratio2
                                        ])
                    
    return overlaps

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description = 'This script extracts reciprocal overlapping structural variants between your study and gnomad ')


    argparser.add_argument('-i', '--input', metavar = 'file', dest = 'bqc_filename', type = str, required = True, help = 'Input VCF file path required')
    argparser.add_argument('-s', '--source', metavar = 'file', dest = 'gnomad_filename', type = str, required = True, help = 'Annotation source VCF file path required')
    argparser.add_argument('-t', '--svtype', metavar = 'name', dest = 'svtype', type = str , required = True, help = 'svtype required') 


    args = argparser.parse_args()

    res = find_overlaps(args.bqc_filename, args.gnomad_filename, args.svtype)
    df = pd.DataFrame(res,  columns = ['SVTYPE', 'BQC_ID', 'CHROM', 'BQC_POS', 'BQC_END', 'BQC_AC', 'BQC_AN', 'BQC_AF','BQC_SVLEN', 'BQC_MISSING', 'gnomad_AC', 'gnomad_AN','gnomad_AF', 'gnomad_SVLEN', 'BQC_overlap', 'gnomad_overlap'] )
    df.to_csv(f"df_{args.svtype}_annotation.tsv", sep="\t", index=False)
    df.to_pickle(f"df_{args.svtype}_annotation.pickle")
