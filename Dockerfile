FROM rocker/r-ver:4.2.2

# This pipeline has no environment definition of its own -- bin/*.R just assumes sceptre, ondisc,
# and the R packages below are already on PATH. It works wherever that happens to already be true
# (e.g. a pre-built HPC environment) and cannot run in an isolated container -- like a Google Batch
# task on Seqera Platform -- without one. This Dockerfile is that environment.
#
# rocker/r-ver ships Debian's reference BLAS/LAPACK, which is dramatically slower than an optimized
# BLAS for the GLM-fitting and resampling work sceptre does per pair -- switched to OpenBLAS below.

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev zlib1g-dev \
    libopenblas0-pthread \
    git ca-certificates && \
    update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 && \
    update-alternatives --set liblapack.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3 && \
    rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("remotes","data.table","dplyr","ggplot2","arrow","BiocManager"), repos="https://cloud.r-project.org")'
RUN R -e 'BiocManager::install("rhdf5lib", update = FALSE, ask = FALSE)'

# Pinned explicitly, as build args with defaults, rather than left to install_github()'s own
# default of "whatever is on the default branch HEAD *at build time*". Without an explicit `ref`,
# rebuilding this exact Dockerfile on two different days could silently install two different
# sceptre versions -- neither pipeline nor sceptre itself would ever tell you that happened, and a
# sceptre_object.rds produced under one version is not guaranteed to be readable, or to mean the
# same thing statistically, under another (see check_sceptre_api.R in EngreitzLab's
# element-gene-power-analysis for what that failure mode looks like when it's guarded against).
# Override with --build-arg to test a different pin; the defaults are what this Dockerfile,
# including the OpenBLAS switch above, was actually validated against.
ARG ONDISC_REF=58e851a95cfd381b316ac40ed4490e569bd40e2f
ARG SCEPTRE_REF=e7866dd4bc158e4415588c4dde0a69da6ac93166
ENV ONDISC_REF=${ONDISC_REF}
ENV SCEPTRE_REF=${SCEPTRE_REF}

RUN R -e 'remotes::install_github("timothy-barry/ondisc", ref = Sys.getenv("ONDISC_REF"), upgrade = "never")'
RUN R -e 'remotes::install_github("katsevich-lab/sceptre", ref = Sys.getenv("SCEPTRE_REF"), upgrade = "never")'
RUN R -e '
library(sceptre); library(ondisc); library(data.table); library(dplyr); library(ggplot2); library(arrow)
cat("sceptre:", as.character(packageVersion("sceptre")), "\n")
cat("ondisc:  ", as.character(packageVersion("ondisc")), "\n")
'
