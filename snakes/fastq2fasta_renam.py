#!/usr/bin/env python

import sys

# obtener archivo fasta y renombrar secuencias

import sys
import os
import gzip

try:
    fastq_file = sys.argv[1]
    fasta_file = sys.argv[2]
except:
    print("""
Obtener archivo fasta y renombrar secuencias de un archivo fastq.
Ejemplo:
fastq2fasta_renam.py <infile.fastq.gz> <outfile.fasta>
    """)
    sys.exit()

# funciones

def fastq_to_fasta(fastq_file, fasta_file, prefix):
    """Convertir un archivo fastq en fasta y cambiar los encabezados."""
    with gzip.open(fastq_file, 'rb') as fq, open(fasta_file, 'w') as fa:
        sequence_number = 1
        
        while True:
            # leer cuatro lineas del archivo fastq
            header = fq.readline().strip()
            if not header:  # si el encabezado esta vacio, es el fin
                break
            sequence = fq.readline().strip()
            fq.readline()  # saltar la linea '+'
            fq.readline()  # saltar la linea de calidad
            
            # crear un nuevo encabezado
            new_header = f">{prefix}_{sequence_number}"
            fa.write(f"{new_header}\n")
            fa.write(f"{sequence.decode("utf-8")}\n")
            
            sequence_number += 1 

# main
prefix = fastq_file.split('/')[-1].strip('.fastq.gz')
fastq_to_fasta(fastq_file, fasta_file, prefix)
