FROM rocker/r-ver:4.2.2

# Environment for running this pipeline in a container (e.g. Google Batch on Seqera), since it has
# none of its own -- see the PR description for details, including why OpenBLAS is switched in below.

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

# Pinned explicitly as build args (override with --build-arg) rather than left to
# install_github()'s default of tracking main's HEAD at build time -- see the PR description.
ARG ONDISC_REF=58e851a95cfd381b316ac40ed4490e569bd40e2f
ARG SCEPTRE_REF=e7866dd4bc158e4415588c4dde0a69da6ac93166
ENV ONDISC_REF=${ONDISC_REF}
ENV SCEPTRE_REF=${SCEPTRE_REF}

RUN R -e 'remotes::install_github("timothy-barry/ondisc", ref = Sys.getenv("ONDISC_REF"), upgrade = "never")'
RUN R -e 'remotes::install_github("katsevich-lab/sceptre", ref = Sys.getenv("SCEPTRE_REF"), upgrade = "never")'
RUN R -e 'library(sceptre); library(ondisc); library(data.table); library(dplyr); library(ggplot2); library(arrow); cat("sceptre:", as.character(packageVersion("sceptre")), "\nondisc: ", as.character(packageVersion("ondisc")), "\n")'
