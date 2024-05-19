# FASTER-REPORT
## Generate reports about fastq files 

This repository contains code for generating html reports for fastq files. It uses the 
`faster` and `faster2` programs to analyise data and a Rmarkdown template to generate the report.
All major platforms are supported - Illumina, PacBio, ONT. The reports are also different depending on the platform selected.

## Running
Two of the many ways to run the code:

- git clone the repo, use conda to create an environment using the `environment.yml` file. After that, execute `faster-report.R`
- as above, but open `faster-report.Rmd` in Rstudio and select *Knit with Parameters*
- use docker, with the helper script `faster-report-docker.sh`

Example run:
```

```