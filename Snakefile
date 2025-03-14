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
        expand(str(outdir) + "/denoised/{lbase}.derep.fna", lbase=LOCBASE),
        expand(str(outdir) + "/denoised/{lbase}.deno.fna", lbase=LOCBASE),
        expand(str(outdir) + "/denoised/{lbase}.nonchimeras.fna", lbase=LOCBASE),
        str(outdir) + "/swarm/coiner.fna",
        str(outdir) + "/swarm/coiner.derep.fna",
        str(outdir) + "/swarm/coiner.swarm13.fna",
        expand(str(outdir) + "/annot/ac_slice_{slice_num}.btout", slice_num=range(1, 11)),
        str(outdir) + "/annot/coiner.btout",
        str(outdir) + "/coiner_otu.tsv",
        str(outdir) + "/coiner_tax.tsv"


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
        
rule dereplicate:
   # start the clustering by removing redundancy
    conda:
        "snakes/vsearch.yaml"
    input:
        str(outdir) + "/merged/{lbase}.fna"
    output:
        str(outdir) + "/denoised/{lbase}.derep.fna",
        str(outdir) + "/denoised/{lbase}.derep.uc"
    threads: cthreads
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --sizeout --uc {output[1]} --threads {threads}
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
    threads: cthreads
    shell:
        """
        vsearch --derep_fulllength {input} --output {output[0]} --sizein --sizeout --uc {output[1]} --threads {threads}
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

rule get_annot_db:
   # download coi reference database
    output:
        str(outdir) + "/annot/MIDORI2_LONGEST_NUC_SP_GB263_CO1_BLAST.fasta.zip"
    params:
        str(outdir) +'/annot'
    shell:
        """
        wget --directory-prefix={params} https://www.reference-midori.info/download/Databases/GenBank263_2024-10-13/BLAST_sp/longest/fasta/MIDORI2_LONGEST_NUC_SP_GB263_CO1_BLAST.fasta.zip
        """

rule unzip_annot_db:
   # download coi reference database
    input:
        str(outdir) + "/annot/MIDORI2_LONGEST_NUC_SP_GB263_CO1_BLAST.fasta.zip"
    output:
        str(outdir) + "/annot/MIDORI2_LONGEST_NUC_SP_GB263_CO1_BLAST.fasta"
    params:
        str(outdir) +'/annot'
    shell:
        """
        unzip {params}/MIDORI2_LONGEST_NUC_SP_GB263_CO1_BLAST.fasta.zip -d {params}
        """

rule build_bl_db:
   # build bast db
    conda:
        "snakes/blast.yaml"
    input:
        str(outdir) + "/annot/MIDORI2_LONGEST_NUC_SP_GB263_CO1_BLAST.fasta"
    output:
        str(outdir) + "/annot/midori.ndb"
    params:
        str(outdir) +'/annot'
    shell:
        """
        makeblastdb -in {input} -dbtype 'nucl' -out {params}/midori
        """

rule split_fasta:
    # split query fasta to speed up the blast search
    input:
        str(outdir) + "/swarm/coiner.swarm13.fna"
    output:
        expand(str(outdir) + "/annot/ac_slice_{slice_num}.fasta", slice_num=range(1, 11))
    params:
        str(outdir) + "/annot"
    shell:
        """
        python3 snakes/slice_fasta.py {input} {params}/ac
        """

rule blast:
    # run blast searches
    conda:
        "snakes/blast.yaml"
    input:
        str(outdir) + "/annot/ac_slice_{slice_num}.fasta",
        str(outdir) + "/annot/midori.ndb"
    output:
        str(outdir) + "/annot/ac_slice_{slice_num}.btout"
    params:
        db=str(outdir) + "/annot/midori"
    threads: cthreads
    shell:
        """
        blastn -db {params.db} -query {input[0]} -outfmt 6 -max_target_seqs 5 -evalue 1e-5 -num_threads {threads} -out {output}
        """

rule merge_blast:
    # merge blast output into single table
    input:
        expand(str(outdir) + "/annot/ac_slice_{slice_num}.btout", slice_num=range(1, 11))
    output:
        str(outdir) + "/annot/coiner.btout"
    params:
       str(outdir) + "/annot"
    shell:
        """
        cat {params}/ac_slice_*.btout > {output}
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


rule build_tax_table:
    # build tax table based on blast output, seed lengths and cutoffs
    conda:
        "snakes/coiner_pylibs.yaml"
    input:
        str(outdir) + "/annot/coiner.btout",
        str(outdir) + "/swarm/coiner.swarm13.fna",
    output:
        str(outdir) + "/coiner_tax.tsv"
    params:
       str(outdir) + "/coiner"
    shell:
        """
        python snakes/build_tax_table.py {input[0]} {input[1]} {params}
        """

# END
