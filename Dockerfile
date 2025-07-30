# To build this image run the following commands:
# docker build --no-cache -t sarand-dev .

# To run the image:
# docker run -v /tmp/sarand:/tmp/sarand -v -it sarand-dev bash

FROM mambaorg/micromamba:1.4-bullseye-slim
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG VER

USER root

# Set the conda environment names to run dependencies
ENV CONDA_BAKTA_NAME='bakta'
ENV CONDA_GRAPH_ALIGNER_NAME='graphaligner'
ENV CONDA_BANDAGE_NAME='bandage'
ENV CONDA_RGI_NAME='rgi'
ENV CONDA_EXE_NAME='micromamba'

# Configure the bakta database
ENV BAKTA_DB='/bakta/db-light'

# Install wget and create the conda environments
RUN apt-get update -y -m && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests -y \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    micromamba create -n ${CONDA_RGI_NAME} -c conda-forge -c bioconda -c defaults -y  \
        rgi=5.2.0 && \
    micromamba create -n ${CONDA_GRAPH_ALIGNER_NAME} -c conda-forge -c bioconda -c defaults -y  \
        graphaligner=1.0.17 && \
    micromamba create -n ${CONDA_BANDAGE_NAME} -c conda-forge -c bioconda -c defaults -y  \
            bandage=0.8.1 && \
    micromamba create -n ${CONDA_BAKTA_NAME} -c conda-forge -c bioconda -c defaults -y \
        bakta=1.8.1 && \
    micromamba create -n sarand -c conda-forge -c bioconda -c defaults -y \
        blast=2.14.0 \
        dna_features_viewer=3.1.2 \
        numpy \
        matplotlib-base \
        gfapy=1.2.3 \
        cd-hit=4.6.8 \
        networkx \
        gzip \
        pandas \
        python \
        pillow \
        biopython

# Install the bakta database
RUN mkdir -p /bakta
WORKDIR /bakta
RUN wget https://zenodo.org/record/7669534/files/db-light.tar.gz?download=1 &&  \
    tar -xzvf db-light.tar.gz?download=1 && \
    rm db-light.tar.gz?download=1 && \
    micromamba run -n ${CONDA_BAKTA_NAME} amrfinder_update --force_update --database ${BAKTA_DB}/amrfinderplus-db

# Copy across the sarand files and install
# Comment this out and replace it with the commented section below once in PyPI
RUN mkdir -p /tmp/sarand
COPY ./sarand /tmp/sarand/sarand
COPY ./setup.py /tmp/sarand/setup.py
COPY ./README.md /tmp/sarand/README.md
RUN micromamba run -n sarand pip install /tmp/sarand && \
    micromamba clean --all --yes && \
    rm -rf /tmp/sarand

# Uncomment when sarand is in PyPI
#RUN micromamba run -n sarand pip install sarand==${VER} && \
#    micromamba clean --all --yes

# Prefix all docker commands with sarand, the user will only need to start from -i
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "micromamba", "run", "-n", "sarand", "sarand"]
CMD ["-h"]
