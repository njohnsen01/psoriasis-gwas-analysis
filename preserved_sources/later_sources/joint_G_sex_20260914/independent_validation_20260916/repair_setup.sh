#!/bin/bash
set -e
D="$(dirname "$0")"
python3 -m pip install --target "$D/python_tools" cmake
export PYTHONPATH="$D/python_tools"
export PATH="$D/python_tools/bin:$PATH"
Rscript -e 'a<-commandArgs(TRUE);install.packages("logistf",lib=a[1],repos="https://cloud.r-project.org",Ncpus=2);.libPaths(c(a[1],.libPaths()));stopifnot(requireNamespace("logistf",quietly=TRUE));cat("LOGISTF_READY",as.character(packageVersion("logistf")),"\n")' "$D/Rlib"
