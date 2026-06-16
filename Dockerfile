# To build this image run the following commands:
# docker build --no-cache -t sarand-dev .

# To run the image:
# docker run -v /tmp/sarand:/tmp/sarand -v -it sarand-dev bash

FROM mambaorg/micromamba:1.4-bullseye-slim
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG VER

USER root

# Create a single conda environment containing all of sarand's dependencies
RUN micromamba create -n sarand -c conda-forge -c bioconda -c defaults -y \
        blast=2.14.0 \
        bandage=0.8.1 \
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
        biopython \
        pyrodigal

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
