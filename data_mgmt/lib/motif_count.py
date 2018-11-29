#!/usr/bin/python

import os
import sys
import textwrap
import argparse
import itertools
import timeit
import time
import csv
import pandas as pd
# import numpy as np
from subprocess import call
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from pathlib import Path


###############################################################################
# Parse arguments
###############################################################################
# start = timeit.default_timer()

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="input fasta file",
                    required=True,
                    type=str,
                    default=sys.stdin)

parser.add_argument("-m", "--motifs",
                    help="input motif list",
                    required=True,
                    type=str)

parser.add_argument("-o", "--output",
                    help="output directory",
                    required=True,
                    type=str)

parser.add_argument("-c", "--chromosome",
                    help = "chromosome",
                    required = False,
                    type = int,
                    default = 0)

parser.add_argument("-b", "--bins",
                    help = "BED file defining the bins",
                    required = False,
                    type = str,
                    default = None)

args = parser.parse_args()

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

# CHANGE THIS TO PASS IN A FILE OBJECT!
def out_results(outfile, outDict, bin_num = None):
    writer = csv.writer(outfile, delimiter = '\t')
    if bin_num:
        for key, value in outDict.items():
            writer.writerow([key] + [value] + [bin_num])
    else:
        for key, value in outDict.items():
            writer.writerow([key] + [value] + [bin_num])

def motif_counter(motif_dict, seq, adj):
    mer = (adj * 2) + 1
    for i in range(0, (len(seq) - mer + 1)):
        motif = seq[i:(i + mer)]
        if motif in motif_dict: motif_dict[motif] += 1

motif_dict = {}
with open(args.motifs) as f:
    motif_list = f.read().splitlines()
adj = int((len(motif_list[0]) - 1) / 2)

if args.bins:
    print("Counting motifs in bins")
    bins = pd.read_table(args.bins, names = ['Chr', 'Start', 'End'])
    if args.chromosome == 0:
        sys.exit("Need to provide a chromosome number along with bins")
    bins = bins[bins['Chr']==('chr'+str(args.chromosome))]
    bins['BIN'] = (bins['End']/1000000).astype(int)
    outName = 'chr' + str(args.chromosome) + '_bins.csv'
    outfile_name = Path(args.output) / outName
    outfile = open(outfile_name, 'w')
    fasta_reader = Fasta(args.input, read_ahead=10000)
    seq = fasta_reader[str(args.chromosome)]
    seqstr = seq[0:len(seq)].seq
    print("Sequence loaded!")
    for index, row in bins.iterrows():
        print("Counting in bin " + str(row['BIN']))
        for m in motif_list:
            motif_dict[m] = 0
        start = row['Start']
        if start > 0: start -= adj
        end = row['End'] + 1 
        if end < len(seq): end += adj
        motif_counter(motif_dict, seqstr[start:end], adj)
        out_results(outfile, motif_dict, row['BIN'])
    outfile.close()
else:
    for m in motif_list:
        motif_dict[m] = 0
    fasta_reader = Fasta(args.input, read_ahead=10000)
    count = 0
    for key in fasta_reader.keys():
        seq = fasta_reader[key]
        seqstr = seq[0:len(seq)].seq

        # print("counting subtypes in record")
        for m in motif_dict.keys():
            occ = occurrences(seqstr, m)
            motif_dict[m] += occ
        count += 1
    outName = 'full_bins.csv'
    outfile_name = Path(args.output) / outName
    outfile = open(outfile_name, 'w')
    writer = csv.writer(outfile, delimiter = '\t')
    for key, value in motif_dict.items():
        writer.writerow([key] + [value])
