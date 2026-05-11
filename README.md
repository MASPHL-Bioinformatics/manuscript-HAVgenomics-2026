This repository contains the analysis code for "Genomic epidemiology as a tool for understanding drivers of hepatitis A community outbreaks in Massachusetts and New Hampshire" by Krasilnikova L, Bouton L, Brock-Fisher TM, Decker E, Godec M, and Thompson Z, et el. All code is organised as [Quarto](https://quarto.org) (`.qmd`) files written mostly in R and one in Python, with a few supporting R helper scripts.

## Getting Started
To run the code, you will need system R and Python installed. Notebooks have been tested with R 4.5.0 and 4.5.2 and various Python versions from 3.11.9-3.12.8.

Install Quarto following the instructions [here](https://quarto.org/docs/get-started/), then clone this repository.

Next, run the wrapper script from within the project directory.

**NOTE**: the wrapper script installs system dependencies and requires sudo access. 
If you don't have sudo privileges, dependencies can be manually installed prior to  running `quarto render` from within the project directory.

```bash
bash source_run.sh
```
The wrapper script will check for Quarto, Python, and R and exit if it can't find all three.
It installs/upgrades necessary system dependencies - this can take a while.
It then checks to make sure all R and Python dependencies are installed, and if not, `sudo` installs them. 
It also uses Quarto to install Tinytex.
The wrapper script renders the Quarto project file-by-file.
In addition to project output files, the wrapper script writes package version information to session info files.

Output files be written to the `outputs` subdirectory within the project folder. This is specified in both `_quarto.yml` and the `_environment` file.

There may still be some Windows line endings (CRLF) around - you can remove them from any given file with `sed -i 's/\r//' FILENAME`

## Project Structure

| File | Description |
|------|-------------|
| `source_run.sh` | Wrapper script — renders the project and collects session info |
| `_environment` | Defines input file paths and environment variables, sourced and exported to all child `.qmd` runs |
| `_quarto.yml` | Project configuration — output directory, execution behaviour, and optional per-file render control |
| `*.qmd` | Quarto markdown files combining narrative and executable code chunks; most execute code without producing a rendered output file |
| `*.R` | Helper scripts sourced by one or more `.qmd` files |

Please note that the JUNIPER helper script sources a C++ file, `cpp_subroutines.cpp`. You will need to install the `Rcpp` package in R: 
`install.packages("Rcpp")`
You also need to install a C++ compiler on your system, such as the version of Rtools specific to your R installation [CRAN RTools](https://cran.r-project.org/bin/windows/Rtools/).

## Rendering Specific Files

By default, the project will render all *.qmd and associated scripts within this project directory. If you want to test or render only a specific subset of files, list them in the `render` queue in `_quarto.yml`, under the project header like so:

```yaml
project:
  render:
    - fig1.qmd
    - fig2.qmd
    - tableS1.R
```

Since each figure is a self-contained qmd file, you can also use a qmd file like a notebook and hardcode in your input paths.

## Notes

Input files are provided in the manuscript to the extent allowed by the IRB. Please contact the authors with any questions.

Delta statistic information: Borges et al 2018. Measuring phylogenetic signal between categorical traits and phylogenies. [DOI LINK](https://doi.org/10.1093/bioinformatics/bty800). `code.R` is a direct copy from [mrborges23/delta_statistic](https://github.com/mrborges23/delta_statistic) and included here for convenience only.

JUNIPER preprint: Specht et al 2025. JUNIPER: Reconstructing Transmission Events from Next-Generation Sequencing Data at Scale. [DOI LINK](https://doi.org/10.1101/2025.03.02.25323192)
