FROM --platform=linux/amd64 mattocci/cmdstanr

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    libglpk-dev

# Install BiocManager and remotes for installing BioC + GitHub packages
RUN Rscript -e "install.packages(c('BiocManager', 'remotes'))"

# Install Bioconductor packages first
RUN Rscript -e "BiocManager::install(c('GenomicRanges', 'IRanges', 'Biobase', 'Rhtslib','Rsamtools'))"

# Install CRAN packages
RUN Rscript -e "install.packages(c( \
    'bayesplot', 'cli', 'dplyr', 'ggh4x', 'ggplot2', 'ggrepel', 'gridExtra', \
    'loo', 'magrittr', 'patchwork', 'ppclust', 'rlang', 'TailRank','scales', 'stringr', \
    'tibble', 'tidyr' \
))"


# Install GitHub packages (cmdstanr, CNAqc, mobster)
RUN Rscript -e "remotes::install_github('caravagn/CNAqc')"
RUN Rscript -e "remotes::install_github('caravagnalab/ctree')"
RUN Rscript -e "remotes::install_github('caravagn/mobster')"
RUN Rscript -e "remotes::install_github('caravagnalab/tickTack', dependencies = FALSE, upgrade = 'never')"

ENV CMDSTAN=/opt/cmdstan/cmdstan-2.29.0
RUN echo "cmdstanr::set_cmdstan_path(Sys.getenv('CMDSTAN'))" >> /root/.Rprofile
# Built-in base packages (already included with R): stats, utils

