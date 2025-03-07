#!/usr/bin/env python

import sys
import pandas as pd

try:
    derepfilename = sys.argv[1]
    swarmfilename = sys.argv[2]
    outfilepref = sys.argv[3]
    
except:
    print("""
Construir tabla de conteos de OTUs a partir de dereplicaciones y swarm.
Ejemplo:
build_otu_table.py <dereplication.uc> <swarm.txt> <oufile_prefix>
    """)
    sys.exit()


# FUNCTINONS

def get_ins_size_derepclusters(derep):
    ins_size_d = {} # diccionario de tamanios en cada muestra
    derepclusters_d = {}
    for i in derep.loc[derep['record_type']!='C'].index.values: # un a fila por secuencia
        crtype = derep['record_type'][i] # tipo de registro actual
        cseqlabel = derep['label_query'][i] # encabezado de secuencia actual
        csize = int(cseqlabel.split('=')[1]) # tamanio en su respecitva muestra
        crawid = cseqlabel.split(';')[0] # id crudo actual
        ccentrawid = derep['label_centroid'][i].split(';')[0] # id crudo del centroide
        
        ins_size_d[crawid] = csize # cada id crudo apunta a su abundancia por muestra
        if crtype == 'S':
            derepclusters_d[crawid] = [crawid] # si es centroide, crea una entrada
        else:
            derepclusters_d[ccentrawid].append(crawid) # de otra forma agrega el id
            
    return ins_size_d, derepclusters_d
        
def get_swarms(swarmfilename, derepclusters_d):
    swarms = []
    with open(swarmfilename, 'r') as infile:
        for line in infile:
            line = line.strip()
            linesep = line.split(' ')
            swarms.append([cclust.split(';')[0] for cclust in linesep]) # guardar solo ids crudos

    # cada lista es un cluster, por ahora cada elemento es un centroide de la dereplicacion
    fullswarms = [] # clusters incluyendo cada centroide por muestra
    for cswarm in swarms:
        cfulls = [] # swarm completo actual
        for ccentroid in cswarm: # por cada centroide actual
            cfulls += derepclusters_d[ccentroid]
        fullswarms.append(cfulls) # agrega el cluster actual completo
        
    return fullswarms

def get_sample_names(derep):
    allseqnames = list(derep['label_query'].unique())
    samplenames = list(set([cname.split('_')[0] for cname in allseqnames])) # suponemos que el formato es: muesta_#deseq;...
    return samplenames
    

# MAIN

derep=pd.read_csv(derepfilename, sep='\t', header=None)
derep.columns=['record_type', 'cluster_number', 'length_size', 'psimil', 'orientation',
               '6', '7', '8', 'label_query', 'label_centroid']

ins_size_d, derepclusters_d = get_ins_size_derepclusters(derep)

fullswarms = get_swarms(swarmfilename, derepclusters_d)

sample_names = get_sample_names(derep)

bopic_otu = pd.DataFrame(0, index=sample_names, columns = [cswarm[0] for cswarm in fullswarms])
for cfswarm in fullswarms:
    cseed = cfswarm[0] # semilla actual = indice de la fila actual
    for cmember in cfswarm:
        samplename = cmember.split('_')[0]  # suponemos que el formato es: muesta_#deseq;...
        bopic_otu.loc[samplename, cseed] += ins_size_d[cmember]

bopic_otu.to_csv(outfilepref +'_otu.tsv', sep='\t')
