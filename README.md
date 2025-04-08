# coiner

Metabarcoding with cox1 amplicon sequences.

# installation

Pull the repository:

`git clone https://github.com/miferg/coiner.git`

Set up and activate a conda environment with snakemake:

`conda create -c conda-forge -c bioconda -n snakemake snakemake`

All dependencies will be installed with conda when the pipeline runs for the first time.

# usage

Store all your fastq files in a same directory. Files must be paired fastq files.

Example:

`snakemake --cores 4 --use-conda --config querydir="my_reads" outdir="coiner_out" splitkey="__COI_R" threads_per_job=2`


