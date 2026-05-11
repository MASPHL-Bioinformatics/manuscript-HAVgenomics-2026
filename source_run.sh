#!/bin/bash
set -e # errexit, exit if a command returns non-zero exit status

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()    { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

# ── Environment ───────────────────────────────────────────────────────────────
[[ -f _environment ]] || error "_environment file not found."
set -a; source _environment; set +a # allexport, automatic export of all defined variables/functions to the environment so subshells/children can access them

# ── Check for _quarto.yml ─────────────────────────────────────────────────────
[[ -f _quarto.yml ]] || error "_quarto.yml not found in current directory."

# ── Output directory ──────────────────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR:?OUTPUT_DIR is not set in _environment}"

# ── Check Quarto install ──────────────────────────────────────────────────────
if command -v quarto &>/dev/null; then
    QUARTO_VERSION=$(quarto --version)
    info "Quarto found: version ${QUARTO_VERSION}"
else
    error "No Quarto installation found. Install from https://quarto.org/docs/get-started/"
fi

# ── Check for R ───────────────────────────────────────────────────────────────
if command -v Rscript &>/dev/null; then
    R_VERSION=$(Rscript -e 'cat(as.character(getRversion()))')
    info "R found: version ${R_VERSION}"
else
    warn "R not found on PATH. R-based documents will fail to render."
    warn "Install R from https://cran.r-project.org/"
    error "R is required. Aborting."
fi

# ── Check for Python ──────────────────────────────────────────────────────────
PYTHON_CMD=""
for cmd in python3 python; do
    if command -v "$cmd" &>/dev/null; then
        PYTHON_CMD="$cmd"
        PYTHON_VERSION=$("$cmd" --version 2>&1)
        info "Python found: ${PYTHON_VERSION}"
        break
    fi
done

if [[ -z "$PYTHON_CMD" ]]; then
    warn "Python not found on PATH. Python-based documents will fail to render."
    warn "Install Python from https://www.python.org/downloads/"
fi

# ── Install system dependencies ───────────────────────────────────────────────
info "Installing system dependencies..."
if command -v apt-get &>/dev/null; then
    sudo apt-get update -qq
    sudo apt-get install -y \
        libfontconfig1-dev \
        libcairo2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libudunits2-dev \
        cmake
    info "System dependencies installed."
else
    warn "apt-get not found — skipping system dependency install. If R package compilation fails, install system libraries manually."
fi

# ── Install R packages ────────────────────────────────────────────────────────
info "Checking R packages..."

sudo Rscript -e '
  ## Install cli first — other packages depend on it and a corrupt install breaks everything
  install.packages("cli", repos = "https://cloud.r-project.org")
  cran_packages <- c(
    "ape", "cli", "coda", "cowplot", "data.table",
    "dplyr", "expm", "extraDistr", "finalfit", "geodist", "ggpattern", "ggpmisc", "ggplot2",
    "ggpubr", "ggraph", "gridExtra", "igraph", "kableExtra", "knitr", "lmPerm",
    "phangorn", "phytools", "plyr", "readxl", "rmarkdown", "sf",
    "splitstackshape", "stringi", "stringr", "tibble", "tidyr",
    "tidyverse", "tinytex", "vegan", "VennDiagram"
  )
  missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  } else {
    message("All CRAN packages already installed.")
  }

  # Always reinstall BiocManager to avoid stale versions, then migrate to 3.22
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(version = "3.22", ask = FALSE)
  bioc_packages <- c("ggtree", "ggtreeExtra", "finalfit")
  missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, version = "3.22")
  } else {
    message("All Bioconductor packages already installed.")
  }

  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")

  github_packages <- list(
    TransPhylo = "xavierdidelot/TransPhylo",
    juniper0   = "broadinstitute/juniper0"
  )
  for (pkg in names(github_packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing GitHub package: ", pkg)
      remotes::install_github(github_packages[[pkg]], dependencies = TRUE)
    }
  }
'
info "R package check complete."

# ── Install Python packages ───────────────────────────────────────────────────
if [[ -n "$PYTHON_CMD" ]]; then
    info "Checking Python packages..."
    # baltic may not be on PyPI — fall back to GitHub install if pip fails
    $PYTHON_CMD -m pip install --quiet pandas matplotlib || \
        warn "Some Python packages failed to install — check pip output above."
    $PYTHON_CMD -m pip install --quiet baltic 2>/dev/null || \
        $PYTHON_CMD -m pip install --quiet \
            git+https://github.com/evogytis/baltic || \
        warn "baltic could not be installed — install manually."
    info "Python package check complete."
fi

# ── Install tinytex ──────────────────────────────────────────────────────────────

quarto install tinytex --no-prompt

# ── Render ────────────────────────────────────────────────────────────────────
info "Starting quarto render..."
quarto render --execute
info "Render complete."

# ── Find source files ─────────────────────────────────────────────────────────
ALL_FILES=$(find . \
    -not -path './.git/*' \
    -not -path "./${OUTPUT_DIR}/*" \
    \( -name "*.qmd" -o -name "*.R" -o -name "*.r" -o -name "*.py" \) \
    | sort)

if [[ -z "$ALL_FILES" ]]; then
    warn "No source files found — skipping session info."
    info "Done."
    exit 0
fi

# ── R session info ────────────────────────────────────────────────────────────
SESSION_INFO_FILE="${OUTPUT_DIR}/session_info.txt"
info "Collecting R session info → ${SESSION_INFO_FILE}"

# R and .qmd files only (exclude .py)
R_FILES=$(echo "$ALL_FILES" | grep -v '\.py$' || true)

if [[ -z "$R_FILES" ]]; then
    warn "No R/.qmd files found — skipping R session info."
else
    Rscript - <<'RSCRIPT' "$SESSION_INFO_FILE" $R_FILES
args      <- commandArgs(trailingOnly = TRUE)
out_file  <- args[1]
src_files <- args[-1]

sink(out_file)

cat("Session info generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("=", 60), "\n\n")

# ── Base session info ──────────────────────────────────────────────────────
cat("BASE SESSION INFO\n")
cat(strrep("-", 60), "\n")
print(sessionInfo())
cat("\n")

# ── Packages detected in source files ─────────────────────────────────────
cat(strrep("=", 60), "\n")
cat("PACKAGES DETECTED IN SOURCE FILES\n")
cat(strrep("-", 60), "\n")

pkg_pattern <- paste0(
    "(?:library|require)\\((['\"]?)([A-Za-z][A-Za-z0-9._]*)\\1\\)",
    "|",
    "([A-Za-z][A-Za-z0-9._]*)::(?:[A-Za-z0-9._]+)"
)

all_pkgs <- character(0)
file_pkgs <- list()

for (f in src_files) {
    if (!file.exists(f)) next
    lines   <- readLines(f, warn = FALSE)
    lines   <- sub("#.*", "", lines)  # strip comments
    matches <- regmatches(lines, gregexpr(pkg_pattern, lines, perl = TRUE))
    matches <- unlist(matches)

    pkgs <- c(
        sub(".*(?:library|require)\\(['\"]?([A-Za-z][A-Za-z0-9._]*)['\"]?\\).*", "\\1",
            grep("library|require", matches, value = TRUE), perl = TRUE),
        sub("([A-Za-z][A-Za-z0-9._]*)::.*", "\\1",
            grep("::", matches, value = TRUE), perl = TRUE)
    )
    pkgs <- unique(pkgs[nchar(pkgs) > 0])
    file_pkgs[[f]] <- pkgs
    all_pkgs <- c(all_pkgs, pkgs)
    cat(sprintf("\n%s\n", f))
    if (length(pkgs) > 0) cat(paste(" ", pkgs, collapse = "\n"), "\n") else cat("  (none detected)\n")
}

all_pkgs <- sort(unique(all_pkgs))

# ── Installed versions for detected packages ───────────────────────────────
cat("\n", strrep("=", 60), "\n", sep = "")
cat("INSTALLED VERSIONS OF DETECTED PACKAGES\n")
cat(strrep("-", 60), "\n")

ip <- as.data.frame(installed.packages()[, c("Package", "Version")],
                    stringsAsFactors = FALSE)

for (pkg in all_pkgs) {
    row <- ip[ip$Package == pkg, ]
    if (nrow(row) > 0) {
        cat(sprintf("  %-30s %s\n", pkg, row$Version))
    } else {
        cat(sprintf("  %-30s (not installed)\n", pkg))
    }
}

sink()
RSCRIPT

    info "R session info written to ${SESSION_INFO_FILE}"
fi

# ── Python session info ───────────────────────────────────────────────────────
if [[ -z "$PYTHON_CMD" ]]; then
    warn "Skipping Python session info — Python not found."
else
    PYTHON_SESSION_FILE="${OUTPUT_DIR}/session_info_python.txt"
    info "Collecting Python session info → ${PYTHON_SESSION_FILE}"

    # Collect .py files and .qmd files that contain python chunks
    PY_FILES=$(echo "$ALL_FILES" | grep '\.py$' || true)
    PY_QMDS=$(grep -rl '```{python' $(echo "$ALL_FILES" | grep '\.qmd$') 2>/dev/null || true)
    ALL_PY_SOURCES=$(echo -e "${PY_FILES}\n${PY_QMDS}" | grep -v '^$' | sort -u || true)

    if [[ -z "$ALL_PY_SOURCES" ]]; then
        warn "No Python source files or Python chunks in .qmd files found — skipping Python session info."
    else
        $PYTHON_CMD - <<PYSCRIPT "$PYTHON_SESSION_FILE" $ALL_PY_SOURCES
import sys, re, importlib.metadata, importlib.util
from datetime import datetime

out_file  = sys.argv[1]
src_files = sys.argv[2:]

with open(out_file, "w") as fh:
    def p(*args, **kwargs):
        print(*args, **kwargs, file=fh)

    p(f"Session info generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    p("=" * 60)
    p()

    # ── Python version info ────────────────────────────────────────────────
    p("BASE SESSION INFO")
    p("-" * 60)
    p(f"Python version : {sys.version}")
    p(f"Platform       : {sys.platform}")
    p(f"Executable     : {sys.executable}")
    p()

    # ── Detect imports in source files ────────────────────────────────────
    p("=" * 60)
    p("PACKAGES DETECTED IN SOURCE FILES")
    p("-" * 60)

    # matches: import X, import X as Y, from X import ..., from X.Y import ...
    import_pattern = re.compile(
        r'^\s*(?:import|from)\s+([A-Za-z][A-Za-z0-9_]*)'
    )

    all_pkgs = set()
    for f in src_files:
        try:
            with open(f) as src:
                lines = src.readlines()
        except FileNotFoundError:
            continue

        pkgs = set()
        for line in lines:
            line = line.split("#")[0]  # strip comments
            m = import_pattern.match(line)
            if m:
                pkgs.add(m.group(1))

        all_pkgs.update(pkgs)
        p(f"\n{f}")
        if pkgs:
            for pkg in sorted(pkgs):
                p(f"  {pkg}")
        else:
            p("  (none detected)")

    # ── Installed versions ─────────────────────────────────────────────────
    p()
    p("=" * 60)
    p("INSTALLED VERSIONS OF DETECTED PACKAGES")
    p("-" * 60)

    # map top-level import name → installed distribution version
    for pkg in sorted(all_pkgs):
        # skip stdlib modules
        if pkg in sys.stdlib_module_names:
            p(f"  {pkg:<30} (stdlib)")
            continue
        try:
            version = importlib.metadata.version(pkg)
            p(f"  {pkg:<30} {version}")
        except importlib.metadata.PackageNotFoundError:
            # some packages have different dist name vs import name (e.g. PIL → Pillow)
            p(f"  {pkg:<30} (not found — may differ from distribution name)")

PYSCRIPT

        info "Python session info written to ${PYTHON_SESSION_FILE}"
    fi
fi

info "Done."