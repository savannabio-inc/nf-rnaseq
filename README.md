<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-rnaseq_logo_dark.png">
    <img alt="nf-core/rnaseq" src="docs/images/nf-core-rnaseq_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/rnaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/rnaseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnaseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/rnaseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rnaseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1400710-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1400710)[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/rnaseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnaseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnaseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/rnaseq** is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.

![nf-core/rnaseq metro map](docs/images/nf-core-rnaseq_metro_map_grey_animated.svg)

> In case the image above is not loading, please have a look at the [static version](docs/images/nf-core-rnaseq_metro_map_grey.png).

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Auto-infer strandedness by subsampling and pseudoalignment ([`fq`](https://github.com/stjude-rust-labs/fq), [`Salmon`](https://combine-lab.github.io/salmon/))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
5. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
6. Removal of genome contaminants ([`BBSplit`](http://seqanswers.com/forums/showthread.php?t=41288))
7. Removal of ribosomal RNA ([`SortMeRNA`](https://github.com/biocore/sortmerna))
8. Choice of multiple alignment and quantification routes:
   1. [`STAR`](https://github.com/alexdobin/STAR) -> [`Salmon`](https://combine-lab.github.io/salmon/)
   2. [`STAR`](https://github.com/alexdobin/STAR) -> [`RSEM`](https://github.com/deweylab/RSEM)
   3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) -> **NO QUANTIFICATION**
9. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
10. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
11. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
12. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
13. Create bigWig coverage files ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
14. Extensive quality control:
    1. [`RSeQC`](http://rseqc.sourceforge.net/)
    2. [`Qualimap`](http://qualimap.bioinfo.cipf.es/)
    3. [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html)
    4. [`Preseq`](http://smithlabresearch.org/software/preseq/)
    5. [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) -> [`Bracken`](https://ccb.jhu.edu/software/bracken/) on unaligned sequences; _optional_
15. Pseudoalignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/) or ['Kallisto'](https://pachterlab.github.io/kallisto/); _optional_)
16. Present QC for raw read, alignment, gene biotype, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

> **Note**
> The SRA download functionality has been removed from the pipeline (`>=3.2`) and ported to an independent workflow called [nf-core/fetchngs](https://nf-co.re/fetchngs). You can provide `--nf_core_pipeline rnaseq` when running nf-core/fetchngs to download and auto-create a samplesheet containing publicly available samples that can be accepted directly as input by this pipeline.

> **Warning**
> Quantification isn't performed if using `--aligner hisat2` due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2 derived genomic alignments. However, you can use this route if you have a preference for the alignment, QC and other types of downstream analysis compatible with the output of HISAT2.

## Quick start

**Requirements:** [Nextflow](https://www.nextflow.io/) `>=24.04.2`, Java, and a container engine (`-profile docker` or `singularity` recommended) or Conda (`-profile conda`).

Validate your install with the bundled test profile (downloads reference data and runs a small public dataset):

```bash
cd nf-rnaseq
nextflow run main.nf -profile test,docker --outdir results_test
```

For day-to-day runs from this repository:

```bash
nextflow run /path/to/nf-rnaseq/main.nf \
  -profile docker \
  --outdir /path/to/results \
  ...
```

> [!NOTE]
> New to Nextflow? See the [nf-core installation guide](https://nf-co.re/docs/usage/installation). Pass parameters on the CLI or via `-params-file`; custom `-c` config files cannot override pipeline parameters ([details](https://nf-co.re/docs/usage/configuration#custom-configuration-files)).

## Usage

### Input options

You must provide **either** a samplesheet **or** direct FASTQ parameters (not both).

#### 1. Samplesheet (multiple samples)

CSV with header `sample,fastq_1,fastq_2,strandedness`. Paths may be local or cloud (`s3://`, `gs://`, etc.). Repeat the same `sample` on multiple rows to merge technical replicates (lanes/runs).

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,/data/rep1_L1_R1.fastq.gz,/data/rep1_L1_R2.fastq.gz,auto
CONTROL_REP1,/data/rep1_L2_R1.fastq.gz,/data/rep1_L2_R2.fastq.gz,auto
TREATMENT_REP1,/data/treat_R1.fastq.gz,/data/treat_R2.fastq.gz,forward
```

Leave `fastq_2` empty for single-end libraries. Set `strandedness` to `forward`, `reverse`, `unstranded`, or `auto` (Salmon subsampling infers strandness when `auto`).

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --fasta genome.fa \
  --gtf genome.gtf \
  --outdir results
```

Or use a pre-built [iGenomes](https://nf-co.re/usage/reference_genomes) reference:

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --genome GRCh38 \
  --outdir results
```

#### 2. Direct FASTQ parameters (single sample)

For one sample, pass FASTQs without creating a CSV. No samplesheet file is written at runtime.

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--fastq_1` | Yes | Read 1 (`.fastq.gz` / `.fq.gz`) |
| `--sample` | Yes | Sample ID (no spaces) |
| `--fastq_2` | No | Read 2 for paired-end; omit for single-end |
| `--input_strandedness` | No | Default `auto`; same as samplesheet `strandedness` |

**Paired-end:**

```bash
nextflow run main.nf \
  -profile docker \
  --fastq_1 /data/sample_R1.fastq.gz \
  --fastq_2 /data/sample_R2.fastq.gz \
  --sample MY_SAMPLE \
  --genome GRCh38 \
  --outdir results
```

**Single-end:**

```bash
nextflow run main.nf \
  -profile docker \
  --fastq_1 /data/sample.fastq.gz \
  --sample MY_SAMPLE \
  --genome GRCh38 \
  --outdir results
```

More detail: [docs/usage.md](docs/usage.md) (samplesheet rules, strandedness, trimming, aligners).

### Reference genome

Provide **one** of:

- `--genome <ID>` — iGenomes bundle (e.g. `GRCh38`, `GRCm39`)
- `--fasta` + `--gtf` (or `--gff`) — custom reference

Optional: `--transcript_fasta`, pre-built indexes (`--star_index`, `--salmon_index`, …), `--save_reference` to reuse indexes on later runs.

### Alignment and quantification routes

| `--aligner` | Quantification | Notes |
|-------------|----------------|-------|
| `star_salmon` (default) | STAR → Salmon | Recommended default |
| `star_rsem` | STAR → RSEM | |
| `hisat2` | None | Alignment + QC only |

Pseudoalignment-only (no genomic aligner): `--skip_alignment --pseudo_aligner salmon` or `kallisto`.

### Profiles

| Profile | Use case |
|---------|----------|
| `docker` | Local/server with Docker |
| `singularity` | HPC / Apptainer |
| `conda` / `mamba` | No containers |
| `test` | CI smoke test (small dataset) |
| `test_full` | Full test dataset |

Combine profiles: `-profile docker,test`.

### Useful flags

| Flag | Effect |
|------|--------|
| `--skip_qc` | Skip most QC processes |
| `--skip_trimming` | Input already trimmed |
| `--skip_alignment` | Pseudoalignment only |
| `--with_umi` | UMI library protocols |
| `--remove_ribo_rna` | SortMeRNA rRNA removal (needs `--ribo_database_manifest`) |
| `-resume` | Continue a failed run |

Full parameter list: [nf-co.re/rnaseq/parameters](https://nf-co.re/rnaseq/parameters) or `nextflow run main.nf --help`.

## Pipeline layout

Entry point: [`main.nf`](main.nf) → [`workflows/rnaseq/main.nf`](workflows/rnaseq/main.nf).

```
main.nf
├── PREPARE_GENOME          # indexes, filtered GTF, auxiliary references
├── RNASEQ                  # QC → trim → align/quantify → QC reports
└── PIPELINE_INITIALISATION / COMPLETION
```

Input is turned into a sample channel via `getInputSamplesList()` (samplesheet or `--fastq_1` / `--sample` / optional `--fastq_2`).

## Modules and subworkflows

Processes live under [`modules/`](modules/) (nf-core + local) and are composed in [`subworkflows/`](subworkflows/).

### nf-core modules (`modules/nf-core/`)

| Stage | Modules |
|-------|---------|
| Input / QC | `cat/fastq`, `fastqc`, `fastp`, `fq/lint`, `fq/subsample`, `umitools/*`, `umicollapse`, `trimgalore`, `bbmap/bbsplit`, `sortmerna` |
| Alignment | `star/align`, `star/genomegenerate`, `hisat2/*`, `salmon/index`, `salmon/quant`, `rsem/*`, `kallisto/*` |
| BAM processing | `samtools/*`, `picard/markduplicates`, `subread/featurecounts`, `stringtie/stringtie` |
| Coverage | `bedtools/genomecov`, `ucsc/bedclip`, `ucsc/bedgraphtobigwig` |
| QC | `rseqc/*`, `qualimap/rnaseq`, `dupradar`, `preseq/lcextrap`, `kraken2/kraken2`, `bracken/bracken` |
| Pseudoalignment | `custom/tx2gene`, `tximeta/tximport`, `summarizedexperiment/summarizedexperiment` |
| Reporting | `multiqc` |
| Utilities | `gunzip`, `gffread`, `custom/getchromsizes`, `custom/catadditionalfasta`, `untar` |

### Local modules (`modules/local/`)

| Module | Role |
|--------|------|
| `gtf_filter` / `gtf2bed` | GTF cleanup and BED12 for QC |
| `preprocess_transcripts_fasta_gencode` | GENCODE transcript FASTA handling |
| `multiqc_custom_biotype` | Biotype breakdown for MultiQC |
| `rsem_merge_counts` | Merge RSEM count matrices |
| `star_align_igenomes` / `star_genomegenerate_igenomes` | STAR with iGenomes layout |

### Subworkflows

| Path | Purpose |
|------|---------|
| `subworkflows/local/prepare_genome` | Build or stage reference indexes |
| `subworkflows/local/align_star` | STAR alignment workflow |
| `subworkflows/local/quantify_rsem` | RSEM quantification |
| `subworkflows/local/utils_nfcore_rnaseq_pipeline` | Init, validation, direct FASTQ input helpers |
| `subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness` | FASTQ QC, trim, rRNA/contaminant removal, auto-strand |
| `subworkflows/nf-core/fastq_align_hisat2` | HISAT2 alignment |
| `subworkflows/nf-core/fastq_fastqc_umitools_trimgalore` | Trim Galore branch |
| `subworkflows/nf-core/fastq_fastqc_umitools_fastp` | fastp branch |
| `subworkflows/nf-core/quantify_pseudo_alignment` | Salmon / Kallisto pseudoalignment |
| `subworkflows/nf-core/bam_*` | Mark duplicates, RSeQC, UMI dedup, stats, bigWig |

Module metadata (inputs/outputs): each `modules/**/meta.yml`.

## Development and testing

Unit tests use [nf-test](https://www.nf-test.com/) in a self-contained venv (Nextflow + nf-test installed locally; no Conda required):

```bash
./scripts/setup-test-venv.sh    # once
./scripts/run-nf-tests.sh       # direct FASTQ input tests by default
```

See [docs/testing.md](docs/testing.md) for full test commands and environment variables.

Upstream docs: [usage](https://nf-co.re/rnaseq/usage) · [output](https://nf-co.re/rnaseq/output) · [parameters](https://nf-co.re/rnaseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnaseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/rnaseq/output).

This pipeline quantifies RNA-sequenced reads relative to genes/transcripts in the genome and normalizes the resulting data. It does not compare the samples statistically in order to assign significance in the form of FDR or P-values. For downstream analyses, the output files from this pipeline can be analysed directly in statistical environments like [R](https://www.r-project.org/), [Julia](https://julialang.org/) or via the [nf-core/differentialabundance](https://github.com/nf-core/differentialabundance/) pipeline.

## Online videos

A short talk about the history, current status and functionality on offer in this pipeline was given by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) on [8th February 2022](https://nf-co.re/events/2022/bytesize-32-nf-core-rnaseq) as part of the nf-core/bytesize series.

You can find numerous talks on the [nf-core events page](https://nf-co.re/events) from various topics including writing pipelines/modules in Nextflow DSL2, using nf-core tooling, running nf-core pipelines as well as more generic content like contributing to Github. Please check them out!

## Credits

These scripts were originally written for use at the [National Genomics Infrastructure](https://ngisweden.scilifelab.se), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

The pipeline was re-written in Nextflow DSL2 and is primarily maintained by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [Seqera Labs, Spain](https://seqera.io/).

The pipeline workflow diagram was initially designed by Sarah Guinchard ([@G-Sarah](https://github.com/G-Sarah)) and James Fellows Yates ([@jfy133](https://github.com/jfy133)), further modifications where made by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) and Maxime Garcia ([@maxulysse](https://github.com/maxulysse)).

Many thanks to other who have helped out along the way too, including (but not limited to):

- [Alex Peltzer](https://github.com/apeltzer)
- [Colin Davenport](https://github.com/colindaven)
- [Denis Moreno](https://github.com/Galithil)
- [Edmund Miller](https://github.com/edmundmiller)
- [Gregor Sturm](https://github.com/grst)
- [Jacki Buros Novik](https://github.com/jburos)
- [Lorena Pantano](https://github.com/lpantano)
- [Matthias Zepper](https://github.com/MatthiasZepper)
- [Maxime Garcia](https://github.com/maxulysse)
- [Olga Botvinnik](https://github.com/olgabot)
- [@orzechoj](https://github.com/orzechoj)
- [Paolo Di Tommaso](https://github.com/pditommaso)
- [Rob Syme](https://github.com/robsyme)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rnaseq` channel](https://nfcore.slack.com/channels/rnaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/rnaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.1400710](https://doi.org/10.5281/zenodo.1400710)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
