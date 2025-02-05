import subprocess
import glob
import os
import pathlib

#modeldir = Path(config["modeldir"]) # dir with models ending with .cm
querydir = Path(config["querydir"]) # dir with assemblies ending with .fna
#LOCBASE = [x.split('__1_Custom_')[0] for x in querydir.iterdir() if x.is_file() and x.suffix in [".gz"]]
LOCBASE = [x.split('__1_Custom_')[0] for x in os.listdir(querydir) if x.endswith('fastq.gz')]
#MODELSBASE = [x.stem for x in modeldir.iterdir() if x.is_file() and x.suffix in [".cm"]]
outdir = "process_out/"

rule all:
    input:
        expand(outdir + "merged/{lbase}.fastq.gz", lbase=LOCBASE),
        expand(outdir + "merged/{lbase}.fna", lbase=LOCBASE),
        expand(outdir + "denoised/{lbase}.derep.fna", lbase=LOCBASE),
        expand(outdir + "denoised/{lbase}.deno.fna", lbase=LOCBASE)
        #expand(outdir + "extracted/{fnaf}_{cmodel}.fna", fnaf=FNABASE, cmodel=MODELSBASE),
        #expand(outdir + "m8/{fnaf}_{cmodel}.m8", fnaf=FNABASE, cmodel=MODELSBASE),
        #outdir + "m8/merged.m8",
        #outdir + "cmsearch_summary.tab",
        #outdir + "cmsearch_summary.tsv"

rule run_fastp:
   # filter, trim and merge pairs
    conda:
        "snakes/fastp.yaml"
    input:
        str(querydir) + "/{lbase}__1_Custom_1.fastq.gz",
        str(querydir) + "/{lbase}__1_Custom_2.fastq.gz"
    output:
        str(outdir) + "filtered/{lbase}__1_Custom_1.f.fastq.gz",
        str(outdir) + "filtered/{lbase}__1_Custom_2.f.fastq.gz",
        str(outdir) + "merged/{lbase}.fastq.gz",
        str(outdir) + "filtered/{lbase}.trim.json",
        str(outdir) + "filtered/{lbase}.trim.html"
    threads: 4
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} --detect_adapter_for_pe --cut_right -W 4 -M 20 --length_required 240 --thread 4 \
        --merge --merged_out {output[2]} --json {output[3]} --html {output[4]}
        """

rule fastq2fasta_renam:
    # transform merged reads to fasta and rename headers
    input:
        str(outdir) + "merged/{lbase}.fastq.gz"
    output:
        str(outdir) + "merged/{lbase}.fna",
    shell:
        """
        python snakes/fastq2fasta_renam.py {input}  {output}
        """
        
rule dereplicate:
   # start the clustering by removing redundancy
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "merged/{lbase}.fna"
    output:
        str(outdir) + "denoised/{lbase}.derep.fna",
        str(outdir) + "denoised/{lbase}.derep.uc"
    threads: 4
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --sizeout --uc {output[1]} --threads 4
        """

rule denoise:
   # denoise dataset (remove sequencing errors based on models of error frequency)
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "denoised/{lbase}.derep.fna"
    output:
        str(outdir) + "denoised/{lbase}.deno.fna",
        str(outdir) + "denoised/{lbase}.deno.uc"
    threads: 4
    shell:
        """
        vsearch --cluster_unoise {input} --centroids {output[0]} --sizein --sizeout --minsize 2 --uc {output[1]} --threads 4
        """

# END