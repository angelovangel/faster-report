# FASTER-REPORT

Generate interactive quality control HTML reports for FASTQ/BAM files. It supports **Illumina**, **PacBio**, and **Nanopore (ONT)** platforms, with tailored visualizations for each.

The reports are generated using the fast Rust utilities [faster](https://github.com/angelovangel/faster) / [faster2](https://github.com/angelovangel/faster2) and [fastkmers](https://github.com/angelovangel/fastkmers).

---

## Execution Options

### 1. Nextflow Pipeline (Recommended)
Orchestrates BAM-to-FASTQ conversion, header parsing, and report generation in a containerized environment. Only `nextflow` and `docker` needed.

```bash
nextflow run main.nf --reads /path/to/fastq_directory/ [options]
```

* **Header Auto-Detection**: Platform (`type`), run date, flowcell ID, and basecaller model are automatically detected from the first FASTQ file header if not explicitly provided.
* **BAM Support**: `.bam` files in the reads directory are automatically converted to `.fastq` using `samtools`.

#### Key Nextflow Parameters
| Parameter | Description | Default |
| :--- | :--- | :--- |
| `--reads` | Path to folder containing reads (required) | `null` |
| `--type` | Sequencer platform (`illumina`, `ont`, or `pacbio`) | Auto-detected (fallback `ont`) |
| `--subsample` | Fraction of reads to subsample for k-mers (`0.1` to `1.0`) | `1.0` |
| `--outdir` | Output directory to save the report | `output` |
| `--outfile` | Name of the output HTML file | `faster-report.html` |
| `--save_raw` | Save raw CSV tables used for plotting (`true`/`false`) | `false` |

---

### 2. Docker Container
Run the tool inside a pre-built container using the provided wrapper script (which auto-mounts your `$HOME` directory):

```bash
./faster-report-docker.sh -p /path/to/fastq_directory/ [options]
```

Or run via Docker CLI directly:
```bash
docker run -it --mount type=bind,src="$HOME",target="$HOME" -w /path/to/workdir aangeloo/faster-report -p /path/to/fastq_directory/
```

---

### 3. Local CLI Script
Run the R script directly on your host machine:

```bash
./faster-report.R -p /path/to/fastq_directory/ [options]
```

#### Local Prerequisites
* **System Utilities**: `pandoc`, `samtools` (for BAM input), and the compiled binaries for `faster`, `faster2`, and `fastkmers` placed in `bin/` or your system `PATH`.
* **R Packages**: `rmarkdown`, `reactable`, `sparkline`, `htmlwidgets`, `dplyr`, `jsonlite`, `optparse`, `funr`, `bslib`, `bsicons`, `scales`.
* **RStudio**: You can also open `faster-report.Rmd` and select **Knit with Parameters**.

---

## Configuration
Modify [nextflow.config](file:///Users/angeloas/code/faster-report/nextflow.config) to customize execution parameters like allocated CPUs, memory limits, or docker bind mounts.
