import subprocess
import glob
import os
import pathlib

querydir = pathlib.Path(config["querydir"])
outdir = pathlib.Path(config["outdir"])
splitkey = config["splitkey"]
cthreads = config["threads_per_job"]
LOCBASE = [x.split(splitkey)[0] for x in os.listdir(querydir) if x.endswith('fastq.gz')]

rule all:
    input:
        expand(str(outdir) + "/filtered/{lbase}.noprim.1.fastq.gz", lbase=LOCBASE),
        expand(str(outdir) + "/filtered/{lbase}.noprim.2.fastq.gz", lbase=LOCBASE),
        expand(str(outdir) + "/merged/{lbase}.fastq.gz", lbase=LOCBASE),
        expand(str(outdir) + "/merged/{lbase}.fna", lbase=LOCBASE),
        expand(str(outdir) + "/merged/{lbase}.oriented.fna", lbase=LOCBASE),
        expand(str(outdir) + "/denoised/{lbase}.derep.fna", lbase=LOCBASE),
        expand(str(outdir) + "/denoised/{lbase}.deno.fna", lbase=LOCBASE),
        expand(str(outdir) + "/denoised/{lbase}.nonchimeras.fna", lbase=LOCBASE),
        str(outdir) + "/swarm/coiner.fna",
        str(outdir) + "/swarm/coiner.derep.fna",
        str(outdir) + "/swarm/coiner.swarm13.fna",
        str(outdir) + "/coiner_otu.tsv",
        str(outdir) + "/coiner.fna",

rule get_annot_db:
   # download coi reference database
    output:
        str(outdir) + "/annot/eKOI_ver1.fasta"
    params:
        str(outdir) +'/annot'
    shell:
        """
        wget --directory-prefix={params} https://github.com/rubenmiguens/eKOI_taxonomy_database/raw/main/eKOI_ver1.fasta
        """

rule run_cutadapt:
   # remove primer sequences
    conda:
        "snakes/cutadapt.yaml"
    input:
        str(querydir) + "/{lbase}"+ splitkey +"1.fastq.gz",
        str(querydir) + "/{lbase}"+ splitkey +"2.fastq.gz"
    output:
        str(outdir) + "/filtered/{lbase}.noprim.1.fastq.gz",
        str(outdir) + "/filtered/{lbase}.noprim.2.fastq.gz",
        str(outdir) + "/filtered/{lbase}.cutadapt.log"
    shell:
        """
        cutadapt -g $(cat snakes/UEA3F.sequence) -G $(cat snakes/HCO2198R.sequence) --discard-untrimmed -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {output[2]}
        """

rule run_fastp:
   # filter, trim and merge pairs
    conda:
        "snakes/fastp.yaml"
    input:
        str(outdir) + "/filtered/{lbase}.noprim.1.fastq.gz",
        str(outdir) + "/filtered/{lbase}.noprim.2.fastq.gz"
    output:
        str(outdir) + "/filtered/{lbase}.qtrim.1.fastq.gz",
        str(outdir) + "/filtered/{lbase}.qtrim.2.fastq.gz",
        str(outdir) + "/merged/{lbase}.fastq.gz",
        str(outdir) + "/filtered/{lbase}.trim.json",
        str(outdir) + "/filtered/{lbase}.trim.html"
    threads: cthreads
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} --cut_right -W 4 -M 20 --length_required 240 --thread {threads} \
        --merge --merged_out {output[2]} --json {output[3]} --html {output[4]}
        """

rule fastq2fasta_renam:
    # transform merged reads to fasta and rename headers
    input:
        str(outdir) + "/merged/{lbase}.fastq.gz"
    output:
        str(outdir) + "/merged/{lbase}.fna",
    shell:
        """
        python snakes/fastq2fasta_renam.py {input}  {output}
        """

rule orient:
   # orient sequences using the ref db
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/merged/{lbase}.fna",
        str(outdir) + "/annot/eKOI_ver1.fasta"
    output:
        str(outdir) + "/merged/{lbase}.oriented.fna",
    threads: cthreads
    shell:
        """
        vsearch --orient {input[0]} --db {input[1]} --fastaout {output} 
        """

rule dereplicate:
   # start the clustering by removing redundancy
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/merged/{lbase}.oriented.fna"
    output:
        str(outdir) + "/denoised/{lbase}.derep.fna",
        str(outdir) + "/denoised/{lbase}.derep.uc"
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --sizeout --uc {output[1]}
        """

rule denoise:
   # denoise dataset (remove sequencing errors based on models of error frequency)
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/denoised/{lbase}.derep.fna"
    output:
        str(outdir) + "/denoised/{lbase}.deno.fna",
        str(outdir) + "/denoised/{lbase}.deno.uc"
    threads: cthreads
    shell:
        """
        vsearch --cluster_unoise {input} --centroids {output[0]} --sizein --sizeout --minsize 2 --uc {output[1]} --threads {threads}
        """

rule uchime:
   # detect and remove chimeras using de novo approach
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/denoised/{lbase}.deno.fna"
    output:
        str(outdir) + "/denoised/{lbase}.nonchimeras.fna"
    threads: cthreads
    shell:
        """
        vsearch --uchime_denovo {input} --nonchimeras {output[0]} --threads {threads}
        """

rule merge_nonchimeras:
   # merge all nonchimeras
    input:
        expand(str(outdir) + "/denoised/{lbase}.nonchimeras.fna", lbase=LOCBASE)
    output:
        str(outdir) + "/swarm/coiner.fna"
    params:
        str(outdir) + "/denoised/"
    shell:
        """
        cat {params}*.nonchimeras.fna > {output}
        """

rule dereplicate_for_swarm:
   # dereplicate nonchimeras to be able to run swarm
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/swarm/coiner.fna"
    output:
        str(outdir) + "/swarm/coiner.derep.fna",
        str(outdir) + "/swarm/coiner.derep.uc"
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --sizein --sizeout --uc {output[1]}
        """

rule swarm:
   # clustering
    conda:
        "snakes/swarm.yaml"
    input:
        str(outdir) + "/swarm/coiner.derep.fna"
    output:
        str(outdir) + "/swarm/coiner.swarm_out.txt",
        str(outdir) + "/swarm/coiner.swarm13.fna"
    threads: cthreads
    shell:
        """
        swarm  -t {threads} -d 13 {input} -o {output[0]} -z --seeds {output[1]}
        """

rule build_otu_table:
    # build otu count table based on swarm output and dereplication
    conda:
        "snakes/coiner_pylibs.yaml"
    input:
        str(outdir) + "/swarm/coiner.derep.uc",
        str(outdir) + "/swarm/coiner.swarm_out.txt",
    output:
        str(outdir) + "/coiner_otu.tsv"
    params:
       str(outdir) + "/coiner"
    shell:
        """
        python snakes/build_otu_table.py {input[0]} {input[1]} {params}
        """

rule get_fna_representatives:
    # get a copy with clean names from cluster representatives
    input:
        str(outdir) + "/swarm/coiner.swarm13.fna",
    output:
        str(outdir) + "/coiner.fna"
    shell:
        """
        awk -F";" '{{print $1}}' {input} > {output}
        """

# END
