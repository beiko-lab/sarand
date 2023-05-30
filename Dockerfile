# To build this image run the following commands:
# docker build --no-cache -t sarand-dev .

# To run the image:
# docker run -v /tmp/sarand:/tmp/sarand -v /Users/aaron/git/sarand/test/aaron/bakta/db-light:/bakta -it sarand-dev bash

FROM mambaorg/micromamba:1.4-bullseye-slim

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Set this to true so sarand is aware and can modify calls accordingly
ENV IS_DOCKER_CONTAINER=1
ENV BAKTA_DB='/bakta/db-light'



# Install wget
USER root
RUN apt-get update -y -m && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests -y \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    micromamba create -n rgi -c conda-forge -c bioconda -y rgi=6.0.2 && \
    micromamba create -n graphaligner -c conda-forge -c bioconda -y graphaligner=1.0.17 && \
    micromamba create -n bakta -c conda-forge -c bioconda -y bakta=1.8.1 && \
    micromamba create -n sarand -c conda-forge -c bioconda -y \
        blast=2.14.0 \
        dna_features_viewer=3.1.2 \
        numpy \
        matplotlib-base \
        gfapy=1.2.3 \
        pandas \
        python \
        biopython

RUN mkdir -p /bakta
WORKDIR /bakta
COPY db-light.tar.gz /bakta/db-light.tar.gz
RUN tar -xzvf db-light.tar.gz && \
    rm db-light.tar.gz && \
    micromamba run -n bakta amrfinder_update --force_update --database /bakta/db-light/amrfinderplus-db

# Copy across the sarand files and install
RUN mkdir -p /tmp/sarand
COPY ./sarand /tmp/sarand/sarand
COPY ./setup.py /tmp/sarand/setup.py
COPY ./README.md /tmp/sarand/README.md

RUN micromamba run -n sarand pip install /tmp/sarand && \
    micromamba clean --all --yes


ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "micromamba", "run", "-n", "sarand", "sarand"]
CMD ["-h"]
