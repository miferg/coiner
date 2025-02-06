import subprocess
import glob
import os
import pathlib

querydir = Path(config["querydir"])
LOCBASE = [x.split('__1_Custom_')[0] for x in os.listdir(querydir) if x.endswith('fastq.gz')]
outdir = "process_out/"

rule all:
    input:
        expand(outdir + "merged/{lbase}.fastq.gz", lbase=LOCBASE),
        expand(outdir + "merged/{lbase}.fna", lbase=LOCBASE),
        expand(outdir + "denoised/{lbase}.derep.fna", lbase=LOCBASE),
        expand(outdir + "denoised/{lbase}.deno.fna", lbase=LOCBASE),
        expand(outdir + "denoised/{lbase}.nonchimeras.fna", lbase=LOCBASE),
        str(outdir) + "swarm/atlascoi.fna",
        str(outdir) + "swarm/atlascoi.derep.fna",
        str(outdir) + "swarm/atlascoi.swarm13.fna"

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

rule uchime:
   # detect and remove chimeras using de novo approach
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "denoised/{lbase}.deno.fna"
    output:
        str(outdir) + "denoised/{lbase}.nonchimeras.fna"
    threads: 4
    shell:
        """
        vsearch --uchime_denovo {input} --nonchimeras {output[0]} --threads 4
        """

rule merge_nonchimeras:
   # merge all nonchimeras
    input:
        expand(outdir + "denoised/{lbase}.nonchimeras.fna", lbase=LOCBASE)
    output:
        str(outdir) + "swarm/atlascoi.fna"
    params:
        str(outdir) + "denoised/"
    shell:
        """
        cat {params}*.nonchimeras.fna > {output}
        """

rule dereplicate_for_swarm:
   # dereplicate nonchimeras to be able to run swarm
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "swarm/atlascoi.fna"
    output:
        str(outdir) + "swarm/atlascoi.derep.fna",
        str(outdir) + "swarm/atlascoi.derep.uc"
    threads: 4
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --sizein --sizeout --uc {output[1]} --threads 4
        """

rule swarm:
   # clustering
    conda:
        "snakes/swarm.yaml"
    input:
        str(outdir) + "swarm/atlascoi.derep.fna"
    output:
        str(outdir) + "swarm/atlascoi.swarm_out.txt",
        str(outdir) + "swarm/atlascoi.swarm13.fna"
    threads: 8
    shell:
        """
        swarm  -t 8 -d 13 {input} -o {output[0]} -z --seeds {output[1]}
        """




# END