#!/usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO

try:
    btoutfilename = sys.argv[1]
    seqfilename = sys.argv[2]
    outfilepref = sys.argv[3]
    
except:
    print("""
Construir tabla de taxonomia de OTUs a partir de la tabla de blast.
Ejemplo:
build_tax_table.py <blast_hits.btout> <sequences.fna> <oufile_prefix>
    """)
    sys.exit()


# FUNCTINONS

def get_sequence_lengths_with_description(fasta_file):
    """Create a dictionary with the length of each sequence"""
    seq_lengths = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths[record.description] = len(record.seq)
    
    return seq_lengths
    

# MAIN

gminpid_d = {'division': 76.66134848484849, # empiric values after analyzing the full Midori2 database, aerage of all minimal pidents per taxa
 'class': 80.44386608595494,
 'order': 83.15110829427071,
 'family': 85.52851374738837,
 'genus': 90.71593556734628,
 'species': 94} # exception, value found in https://doi.org/10.1673/031.012.1601 rounded up

# load blast output

atlasblast = pd.read_csv(btoutfilename, sep='\t', header=None)
atlasblast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# bet best hits based in bitscore

atlasbest = atlasblast.sort_values('bitscore', ascending=False).groupby('qseqid').first().reset_index().copy()

# get sequence lengths, and query coverage per alignment

seedlens_d = get_sequence_lengths_with_description(seqfilename)
atlasbest['qseq_len'] = [seedlens_d[crow['qseqid']] for ci, crow in atlasbest.iterrows()]
atlasbest['qcov'] = atlasbest['length']/atlasbest['qseq_len']

# remove  alignmens with low qcov

atlasbest = atlasbest.loc[atlasbest['qcov']>=0.6]
atlasbest.shape

# add taxonomy based on means of min pident values

tax_d = {}
for crank in ['division', 'class', 'order', 'family', 'genus', 'species']:
    tax_d[crank] = []

for ci, crow in atlasbest.iterrows():
    taxonpath = crow['sseqid'].split(';')[2:]
    ltax = 'Eukaryota' # initialize the last tax variable in case pident is too low
    for crank in tax_d.keys():
        ctax = taxonpath.pop(0)
        if crow['pident'] >= gminpid_d[crank]:
            tax_d[crank].append(ctax)
            ltax = ctax # current tax becomes last tax in case there are no taxonomic assignments left
        else:
            tax_d[crank].append(ltax +'.U.'+ crank) # this is useful in gappa output

for crank in tax_d.keys():
    atlasbest[crank] = tax_d[crank]

# extract taxonomy from blast table

atlasbest['seed'] = atlasbest['qseqid']
tax_tab = atlasbest[['seed', 'division', 'class', 'order', 'family', 'genus', 'species']]

# write

tax_tab.to_csv(outfilepref +'_tax.tsv', sep='\t', index=False)
