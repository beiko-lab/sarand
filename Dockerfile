FROM condaforge/mambaforge
ARG MAMBA_DOCKERFILE_ACTIVATE=1 

# metadata
LABEL base.image="miniconda3"
LABEL version="0.0.1"
LABEL software="sarand"
LABEL description="Tool to extract AMR neighbourhoods from metagenomic graphs"
LABEL website="https://github.com/maguire-lab/sarand"
LABEL documentation="https://github.com/maguire-lab/sarand/blob/master/README.rst"
LABEL license="https://github.com/maguire-lab/sarand/blob/master/LICENSE"
LABEL tags="Genomics"

# maintainer
MAINTAINER Finlay Maguire <finlaymaguire@gmail.com>

#ADD . /
RUN git clone https://github.com/maguire-lab/sarand
WORKDIR sarand

RUN mamba env create -f conda_env.yaml 
RUN conda run -n sarand pip install .
RUN apt-get update 
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y libglu1-mesa-dev build-essential mesa-common-dev libfontconfig1 software-properties-common zip

RUN wget https://github.com/rrwick/Bandage/releases/download/v0.9.0/Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip && unzip Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip && chmod +x Bandage_Ubuntu-x86-64_v0.9.0.AppImage && cp Bandage_Ubuntu-x86-64_v0.9.0.AppImage /opt/conda/envs/sarand/bin/Bandage

# BandageNG altermative
#RUN wget https://github.com/asl/BandageNG/releases/download/v2022.08/BandageNG-b05263d-x86_64.AppImage && chmod +x BandageNG-b05263d-x86_64.AppImage && cp BandageNG-b05263d-x86_64.AppImage /opt/conda/envs/sarand/bin/Bandage
#RUN apt-get update && apt-get install -y libglx0 libopengl0  libegl1

SHELL ["conda", "run", "-n", "sarand", "/bin/bash", "-c"]
ARG APPIMAGE_EXTRACT_AND_RUN=1
RUN prokka --setupdb
RUN /bin/bash test/test.sh 
ENTRYPOINT ["conda", "run", "-n", "sarand"]
