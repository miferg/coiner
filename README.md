# coiner

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-managed-green.svg)](https://docs.conda.io/en/latest/miniconda.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**coiner** is a metabarcoding pipeline designed for the processing of **COI (Cytochrome c oxidase subunit I)** amplicon sequences.

The pipeline processes paired-end FASTQ files through quality filtering, merging, sequence orientation, and denoising (UNOISE3), culminating in OTU clustering via **Swarm**. It automatically handles the download of the **eKOI** reference database for sequence orientation. 

For post-processing, we recommend using the `DARN` (https://github.com/hariszaf/darn) pipeline to discard non-eukaryotic sequences.



---

## 🛠️ Installation

`coiner` relies on Snakemake to manage the workflow and Conda to manage software dependencies.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/miferg/coiner.git
    cd coiner
    ```

2.  **Create the environment:**
    You only need a base environment with Snakemake.

    ```bash
    conda create -c conda-forge -c bioconda -n coiner-env snakemake
    ```

3.  **Activate:**
    ```bash
    conda activate coiner-env
    ```

*Note: The pipeline will automatically create isolated environments for tools like `vsearch`, `cutadapt`, `fastp`, and `swarm` upon the first run.*

---

## 📂 Input Formatting

To ensure the pipeline recognizes your files correctly, please follow this naming convention:

1.  **Format:** Paired-end FASTQ files (`.fastq.gz`).
2.  **Naming:** Files must end with `1.fastq.gz` (forward) and `2.fastq.gz` (reverse).
3.  **Delimiter:** The filename must contain a delimiter (specified by `splitkey`) that separates the sample name from the rest of the filename.

**Example:**
If `splitkey="_"`:
* `SampleA_1.fastq.gz` & `SampleA_2.fastq.gz` -> Sample ID: **SampleA**
* `Site1-Rep1_LibA_1.fastq.gz` -> Sample ID: **Site1-Rep1**

---

## 🚀 Usage

Run the pipeline from the repository root directory.

### Example Command

```bash
snakemake --cores 8 --use-conda --config querydir="raw_reads/" outdir="coiner_results/" splitkey="_" threads_per_job=4
````

### Configuration Parameters

| Parameter | Description |
| :--- | :--- |
| `--cores` | Total number of CPU cores available to Snakemake. |
| `--use-conda` | **Required.** Enables automatic dependency management. |
| `querydir` | Path to the directory containing your raw paired-end FASTQ files. |
| `outdir` | Path where results and intermediate files will be stored. |
| `splitkey` | The character separating the sample name from the rest of the filename (e.g., `_` or `.`). |
| `threads_per_job` | Number of threads assigned to computationally intensive steps (e.g., `fastp`, `vsearch`, `swarm`). |

-----

## 🔬 Pipeline Workflow

The `coiner` pipeline executes the following steps:

1.  **Database Setup:** Downloads the **eKOI** taxonomy database (ver1) for reference.
2.  **Primer Removal (`cutadapt`):** Removes primers **UEA3F** and **HCO2198R**. Untrimmed reads are discarded.
3.  **QC & Merging (`fastp`):**
      * Quality trimming.
      * Merging paired-end reads (min overlap: 20, min length: 240).
4.  **Orientation (`vsearch`):** Orients merged sequences to the correct strand based on the eKOI reference database.
5.  **Denoising (Per Sample):**
      * **Dereplication:** `vsearch --derep_fulllength`.
      * **Denoising:** `vsearch --cluster_unoise` (UNOISE3 algorithm) to correct sequencing errors.
      * **Chimera Removal:** `vsearch --uchime_denovo`.
6.  **Global Clustering (`swarm`):**
      * Non-chimeric sequences from all samples are pooled.
      * Global dereplication.
      * Clustering using **Swarm** with a difference `d=13`.
7.  **OTU Table Generation:** Maps reads back to clusters to generate the final abundance table.

-----

## 📤 Output

The `outdir` will contain several subdirectories. The most important files are in the root of the output directory:

```text
coiner_out/
├── coiner_otu.tsv        # The final OTU abundance table (samples x OTUs)
├── coiner.fna            # Representative sequences for each OTU
├── swarm/                # Intermediate global clustering files
├── denoised/             # Intermediate per-sample denoised sequences
├── merged/               # Merged and oriented FASTA files
└── filtered/             # Trimmed FASTQ files and QC reports (fastp HTML/JSON)
```

-----

## 📄 License

This project is licensed under the MIT License.
