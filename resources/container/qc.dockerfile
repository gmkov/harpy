FROM mambaorg/micromamba
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    falco \ 
    fastp \ 
    multiqc \ 
    r-base \ 
    r-dplyr \ 
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-knitr \ 
    r-magrittr \ 
    r-plotly \ 
    r-rmarkdown \ 
    r-tidyr \ 
    r-viridislite

RUN micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY workflow/scripts/ /home/

RUN chmod +x /home/scripts/* && cp home/scripts/* $CONDA_PREFIX/bin/