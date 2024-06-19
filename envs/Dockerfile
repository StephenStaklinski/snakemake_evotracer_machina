# hash:sha256:09c569e084991c607e602dbd5438fba69be9a1c4dc397fb009d766ae9380373b
FROM registry.codeocean.com/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN mamba install -y --channel=conda-forge --channel=bioconda --channel=r \
        biopython==1.81 \
        ete3==3.1.3 \
        flash==1.2.11 \
        numpy==1.25.1 \
        openjdk==20.0.2 \
        openssl==3.1.4 \
        pandas==2.0.3 \
        r-gtools==3.9.4 \
        r-rgl==1.2.1 \
        r-statebins==1.4.0 \
        r-tidygraph==1.2.3 \
        scipy==1.11.1 \
        umi_tools==1.1.4 \
    && mamba clean -ya

RUN Rscript -e 'remotes::install_version("RColorBrewer", "1.1-3")' \
    && Rscript -e 'remotes::install_version("TreeTools", "1.10.0")' \
    && Rscript -e 'remotes::install_version("UpSetR", "1.4.0")' \
    && Rscript -e 'remotes::install_version("ade4", "1.7-22")' \
    && Rscript -e 'remotes::install_version("adegenet", "2.1.10")' \
    && Rscript -e 'remotes::install_version("adespatial", "0.3-23")' \
    && Rscript -e 'remotes::install_version("adiv", "2.2")' \
    && Rscript -e 'remotes::install_version("ape", "5.7-1")' \
    && Rscript -e 'remotes::install_version("aplot", "0.2.2")' \
    && Rscript -e 'remotes::install_version("benthos", "1.3-8")' \
    && Rscript -e 'remotes::install_version("betapart", "1.6")' \
    && Rscript -e 'remotes::install_version("broom", "1.0.5")' \
    && Rscript -e 'remotes::install_version("circlize", "0.4.15")' \
    && Rscript -e 'remotes::install_version("colorspace", "2.1-0")' \
    && Rscript -e 'remotes::install_version("data.table", "1.14.8")' \
    && Rscript -e 'remotes::install_version("extrafont", "0.19")' \
    && Rscript -e 'remotes::install_version("forcats", "1.0.0")' \
    && Rscript -e 'remotes::install_version("ggnewscale", "0.4.9")' \
    && Rscript -e 'remotes::install_version("ggraph", "2.1.0")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.4")' \
    && Rscript -e 'remotes::install_version("ggstance", "0.3.6")' \
    && Rscript -e 'remotes::install_version("ggtext", "0.1.2")' \
    && Rscript -e 'remotes::install_version("ggthemes", "4.2.4")' \
    && Rscript -e 'remotes::install_version("ggupset", "0.3.0")' \
    && Rscript -e 'remotes::install_version("graphlayouts", "1.0.2")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("igraph", "1.5.1")' \
    && Rscript -e 'remotes::install_version("lemon", "0.4.7")' \
    && Rscript -e 'remotes::install_version("pacman", "0.5.1")' \
    && Rscript -e 'remotes::install_version("pals", "1.8")' \
    && Rscript -e 'remotes::install_version("phangorn", "2.11.1")' \
    && Rscript -e 'remotes::install_version("pheatmap", "1.0.12")' \
    && Rscript -e 'remotes::install_version("phytools", "2.0-3")' \
    && Rscript -e 'remotes::install_version("quadprog", "1.5-8")' \
    && Rscript -e 'remotes::install_version("readr", "2.1.4")' \
    && Rscript -e 'remotes::install_version("readxl", "1.4.3")' \
    && Rscript -e 'remotes::install_version("reticulate", "1.34.0")' \
    && Rscript -e 'remotes::install_version("rlang", "1.1.2")' \
    && Rscript -e 'remotes::install_version("scales", "1.2.1")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("tibble", "3.2.1")' \
    && Rscript -e 'remotes::install_version("tidygraph", "1.2.3")' \
    && Rscript -e 'remotes::install_version("vegan", "2.6-4")' \
    && Rscript -e 'remotes::install_version("viridis", "0.6.4")'

RUN Rscript -e 'remotes::install_github( \
        "GuangchuangYu/scatterpie", \
        ref = "abecacbc098340d83758c79fb2bacaed19631a6d")' \
    && Rscript -e 'remotes::install_github( \
        "Nowak-Lab/EvoTraceR", \
        ref = "86cf0a1da18eb8803c6716e1c163f38b847f60a9")' \
    && Rscript -e 'remotes::install_github( \
        "hrbrmstr/ggchicklet", \
        ref = "64c468dd0900153be1690dbfc5cfb35710da8183")' \
    && Rscript -e 'remotes::install_github( \
        "igraph/rigraph", \
        ref = "90ca43c92f5c2576059f1a5baa9553c40d508ff3")'

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(c( \
        "Biostrings", \
        "IRanges", \
        "biomaRt", \
        "ggtree", \
        "phyloseq", \
        "treeio" \
    ))' # Original versions: 2.70.1 2.36.0 2.58.0 3.10.0 1.46.0 1.26.0

COPY postInstall /
RUN /postInstall
