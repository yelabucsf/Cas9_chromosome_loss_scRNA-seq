# What’s this repo?

This repository contains all the code to reproduce the results and figures in the paper **Mitigation of chromosome loss in clinical CRISPR-Cas9-engineered T cells, bioRxiv 2023** ([doi.org/10.1101/2023.03.22.533709](https://doi.org/10.1101/2023.03.22.533709)).

## What’s in it?

The repo contains:
1. Scripts for producing the results and figures in the paper:
   - `crop_seq_and_cell_state_analysis.ipynb`
   - `chr14_11gRNA_heatmap.ipynb`
2. Scripts and files for processing the raw fastq and intermediate files to produce the final dataset (gene-expression matrix, sgRNA assignments, copy-number estimates from inferCNV and estimated aneuploidy events):
   - `data_processing/guide_assign_binomial.ipynb`
   - `data_processing/ainfercnv_prep.ipynb`
   - `data_processing/breakpoint_calc.ipynb`
   - `data_processing/inferCNV/Batch{1-9}.R`
   - `data_processing/feature_reference_v2.csv`

Most of the code we provide refers to the main (CROP-seq) screen presented in our paper. The same code applies to the two other screens (Chr14 and CART) with minor changes.

## What’s NOT in it?

The repo doesn’t contain data files. All the raw, intermediate and processed data files are available on the GEO series associated with the project (linked from our publication).

## Required data files

If you only need access to the fully processed data, you will likely want to download the following files from our GEO series:
- `qced.h5ad`
- `concat_InferCNV.pkl`
- `aneuploidy_events.csv`
- `inferCNVgeneName.txt`

# Code for results and figures

Below are the notebooks containing all the code for producing all the results and figures presented in the paper. For each notebook, we provide a high-level description of what the code does, specify which figures are generated by the code, and list all the input data files needed to run the code.

## `crop_seq_and_cell_state_analysis.ipynb`

**Description**: Most of the high-level analysis of the CROP-seq library, including the definition of aneuploidy events (based on the breakpoint estimates), comparison between targeted and lost chromosomes, cell cycle analysis, differential gene expression, and cell state analysis (UMAPs and clusters). Also includes cell state analysis of the TRAC library.

**Produced figures**: Fig. 2C, Fig. 2B, Fig. S3H, data for Fig. 3C (`n_cells_per_chrom_loss_and_cell_cycle_status.csv`), data for Fig. 3A (`dge_results_v3.csv`, `dge_total_gene_scores_v3.csv`), Fig. 3B, Fig. S2C-D, Fig. S4B, Fig. S4A, data for Fig. S2E (`chrom_loss_status_and_cluster_cell_counts.csv`), data for Fig. S4C (`cell_cycle_and_cluster_cell_counts.csv`), Fig. S1C-D, data for Fig. S1E (`chrom_loss_status_and_cluster_cell_counts.csv` - notice that this is a different file than the one used for Fig. S2E)

**Required data files**: for the CROP-seq analysis: `qced.h5ad`, `concat_InferCNV.pkl`, `inferCNVgeneName.txt`, `aneuploidy_events.csv`; for the cell state analysis of the TRAC library: `TRAC_Aneuploidy_events.csv`, `cas9ProcessedAneuploidyStatus.h5ad`

## `chr14_11gRNA_heatmap.ipynb`

**Description**: Performs the inferCNV analysis that is used as the basis for the aneuploidy events for cells treated with the 11 gRNAs targeting the TRAC locus of chromosome 14. 

**Produced figures**: Fig. 1B

**Required data files**: The processed h5ad file, rawSingletForR15Kcells.h5ad and the input file inferCNVgeneName.txt. 

# Re-processing the data files

## Raw data files

The raw data files for the CROP-seq dataset originate from three libraries, each with 24 fastq files:
1. **The gene-expression (GEX) library** contains the files `CTJD02{A-F}_S{1-6}_L00{3/4}_R{1/2}_001.fastq.gz`, where A-F and 1-6 (respectively) correspond to the 6 lanes of the 10x (all on one chip), L003/L004 corresponds to the 2 lanes of the Illumina sequencer, and R1/R2 are the two ends of the paired-end sequencing.
2. **The Multi-Seq library** contains the files `CTJD02{G-L}_S{7-12}_L00{3/4}_R{1/2}_001.fastq.gz`, where G-L and 7-12 (respectively) correspond to the 6 lanes of the 10x (all on one chip), L003/L004 corresponds to the 2 lanes of the Illumina sequencer, and R1/R2 are the two ends of the paired-end sequencing.
3. **The guide/sgRNA library** contains the files `JDCT003{A-F}_S{6-1}_L001_R{1/2}_001.fastq.gz` and `JDCT004{A-F}_S{1-6}_L001_R{1/2}_001.fastq.gz`, where A-F and 6-1 or 1-6 (respectively) correspond to the 6 lanes of the 10x (all on one chip), and R1/R2 are the two ends of the paired-end sequencing. The two JDCT003 and JDCT004 files from the same sample are merged later in the analysis to maximize the number of cells with called guides.

## Step 1: Processing of the GEX and Multi-Seq libraries into cellranger files

We use the `cellranger multi` command to process the GEX and Multi-Seq fastq files:

```
cellranger-7.0.0/bin/cellranger multi --id=sample{1-6} --csv=CTJD02{A-F}{G-L}.csv --localcores=29
```

The command has to be run 6 times, each with a different sample (`--id=sample1`, `--id=sample2`, etc.) and a corresponding CSV file (`--csv=CTJD02AG`, `--csv=CTJD02BH`, etc.).
Each GEX library lane A-F corresponds to Multi-Seq lane G-L (i.e. A goes with G, B goes with H, etc.).

The content of the CSV files `CTJD02{A-F}{G-L}.csv` is:

```
[gene-expression]
reference,/home/ssm-user/references/refdata-gex-GRCh38-2020-A/
cmo-set,featureRefMulti.csv

[libraries]
fastq_id,fastqs,feature_types
CTJD02{A-F},/GEX-fastq-dir/,Gene Expression
CTJD02{G-L},/multiseq–fastq-dir/,Multiplexing Capture

[samples]
sample_id,cmo_ids
sample1,CTJD02A
sample2,CTJD02B
sample3,CTJD02C
sample4,CTJD02D
sample5,CTJD02E
sample6,CTJD02F
```

Here too, make sure each that each of 6 CSV files has the correct `fastq_id` values under `[libraries]` (`CTJD02A` and `CTJD02G`, `CTJD02B` and `CTJD02H`, etc.). The `[samples]` section can remain identical across all 6 files.


All 6 samples use the same multiplexing oligos, provided in `featureRefMulti.csv`. The content for this CSV file is:
```
id,name,read,pattern,sequence,feature_type
CTJD02A,CTJD02A,R2,5P(BC),GGAGAAGA,Multiplexing Capture
CTJD02B,CTJD02B,R2,5P(BC),CCACAATG,Multiplexing Capture
CTJD02C,CTJD02C,R2,5P(BC),TGAGACCT,Multiplexing Capture
CTJD02D,CTJD02D,R2,5P(BC),GCACACGC,Multiplexing Capture
CTJD02E,CTJD02E,R2,5P(BC),AGAGAGAG,Multiplexing Capture
CTJD02F,CTJD02F,R2,5P(BC),TCACAGCA,Multiplexing Capture
```

**Output files**:
For each of the 6 samples, a CSV file named `assignment_confidence_table.csv` is created. In the GEX directory, they are named `sample{1-6}-assignment_confidence_table.csv`.

## Step 2: Processing of the sgRNA and GEX libraries

Here we:
1. Concatenate the fastq files from respective lanes from the two sequencing runs of the sgRNA libraries, specifically JDCT003 and JDCT004. By running the full pipeline on the concatenated guide libraries, we maximize the number of cells with a single, confidently assigned guide.
2. The `cellranger count` command then counts the number of reads per guide per cell and generates the gene-expression matrices.

### 1. Combining guide libraries

The following command concatenates the raw fastq files from the JDCT003 and JDCT004 sequencing runs for each of the 6 samples (A-F) and both paired-ends (R1 and R2):

```bash
cat JDCT003A_S6_L001_R1_001.fastq.gz JDCT004A_S1_L001_R1_001.fastq.gz > JDCT005A_S1_L001_R1_001.fastq.gz
```

(the example is for read R1 of sample A)

### 2. Running cellranger count

The following command is run for each of the 6 samples (you know the drill on how to interpret `{A-F}`):

```bash
cellranger-7.0.1/bin/cellranger count --id=sample{A-F} --transcriptome=/home/ssm-user/references/refdata-gex-GRCh38-2020-A/ --libraries=/home/ssm-user/csvs/sample{A-F}.csv --feature-ref=/data/feature_reference_v2.csv --localcores=62
```

**Input files**:

1. `feature_reference_v2.csv` (provided in `data_processing/` in the repo): Specifies all the ~400 guides used.

2. `sample{A-F}.csv` (example provided for `sampleA.csv`):

    ```csv
    fastqs,sample,library_type,
    /gex_dir/,CTJD02A,Gene Expression,
    /sgrna_dir/,JDCT005A,CRISPR Guide Capture,
    ```

**Output files**:

1. `filtered_feature_bc_matrix.h5`: Gene expression (GEX) matrices for each lane.
    - GEO file names: `sample{A-F}-filtered_feature_bc_matrix.h5`

2. `protospacer_calls_per_cell.csv`: CSV file mapping between cell barcodes and guides, enumerating how many guide counts are found in each cell.
    - Contains the following columns:
        - *cell_barcode*: cell barcode
        - *num_features*: how many different guides were identified
        - *feature_call*: the names of the identified guides
        - *num_umis*: UMIs counted for each guide 
    - GEO file names: `sample{A-F}-protospacer_calls_per_cell.csv`
    
## Step 3: Processing GEX, MultiSeq, and guide calls into a single AnnData object

The notebook `data_processing/guide_assign_binomial.ipynb` integrates GEX, Multi-Seq, and guide calls into one consolidated AnnData object.

**Input files**:

For each of the six samples (`{A-F}`):

1. `filtered_feature_bc_matrix.h5`: Gene expression matrices for each lane from the previous step.

2. `protospacer_calls_per_cell.csv`: A mapping CSV file between cell barcodes and guides from the previous step.

3. `assignment_confidence_table.csv`: Derived from *Step 1*, this table contains data essential for the guide assignment process.

**What the code does**:

- Integrate the raw h5 gene-expression data with the guide counts from the six samples.
- Determine the first and second most common guide per cell.
- Utilize a binomial test to determine whether the most common guide significantly surpasses the counts of the second most common guide.
- Incorporate sample calls into the AnnData object.
- Calculate quality control (QC) metrics.

**Output file**: `fully_processed.h5ad` - a consolidated AnnData file with the processed data (available on GEO).

## Step 4: Constructing gene metadata file for inferCNV

The `gtf_to_position_file.py` script (provided by inferCNV [here](https://github.com/broadinstitute/infercnv/blob/master/scripts/gtf_to_position_file.py)) constructs a metadata file with the chromosomal positions of genes presented in the copy-number matrix of inferCNV.

**Command line**:
```bash
python gtf_to_position_file.py genes.gtf inferCNVgeneName.txt
```

Following the script execution, the ENS gene IDs in `inferCNVgeneName.txt` are translated into gene names.

**Input file**: `genes.gtf` from cellranger's reference (`refdata-gex-GRCh38-2020-A/genes/genes.gtf`)

**Output file**: `inferCNVgeneName.txt` (available on GEO)

## Step 5: Quality control and preparing files for inferCNV

The notebook `data_processing/infercnv_prep.ipynb` applies quality control measures and prepares the files for inferCNV analysis.

**Input file**: `fully_processed.h5ad` - the consolidated AnnData file from *Step 3*.

**What the code does**:
- Implement quality control by applying the singlet filter, a filter for less than 10K total counts, and a filter for less than 10% mitochondrial content.
- To enhance inferCNV's computational efficiency, data is segmented into 9 batches, each representing a different cell subset.

**Output files**:

1. `qced.h5ad`: the final processed h5ad file utilized for most analyses (available on GEO).

2. Files for inferCNV (3 files for each of the 9 batches):
    - `annotations_Batch{1-9}.csv`
    - `genes_Batch{1-9}.csv`
    - `counts_Batch{1-9}.h5ad`
    
## Step 6: Running inferCNV

The R scripts `data_processing/inferCNV/Batch{1-9}.R` conduct the inferCNV analysis for each batch of data.

**Input files**:
1. `annotations_Batch{1-9}.csv` (from *Step 5*).
2. `genes_Batch{1-9}.csv` (from *Step 5*).
3. `counts_Batch{1-9}.h5ad` (from *Step 5*).
4. `inferCNVgeneName.txt` (from *Step 4*).

**What the scripts do**: Process and analyze the segmented data batches using inferCNV to quantify copy number variations along the genome in each cell.

**Output file**: `infercnv.observations.txt` (for each of the 9 batches) (available on GEO as `Batch{1-9}-infercnv.observations.txt`)

## Step 7: Concatenating inferCNV outputs and calculating potential breakpoints

The notebook `data_processing/breakpoint_calc.ipynb` combines the outputs from the inferCNV analysis. The code then scans each chromosome for potential "breakpoints" where the difference between the average inferCNV values to the left and right is greatest.

**Input files**:
1. `qced.h5ad` (from *Step 5*).
2. `infercnv.observations.txt` for each of the 9 batches (from *Step 6*).
3. `inferCNVgeneName.txt` (from *Step 4*).

**What the code does**:
- Read and consolidate the `infercnv.observations.txt` files from the 9 batches.
- Save the combined data as a pickled object for faster future loading.
- Identify and quantify potential chromosomal breakpoints based on the inferCNV values.

**Output files**:
1. `concat_InferCNV.pkl` - the concatenated inferCNV data stored as a pickle object for swift access (available on GEO)
2. `aneuploidy_events.csv` - for each cell and chromosome, the number of inferCNV genes and average inferCNV values to the left and right of the most likely breakpoint (available on GEO).

