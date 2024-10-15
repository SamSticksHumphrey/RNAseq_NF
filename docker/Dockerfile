# Dockerfile for RNAseq Container

FROM condaforge/miniforge3

# basic configuration
RUN apt-get update \
        && apt-get install -y \
                ed \
                less \
                locales \
                gzip \
                unzip \
		pigz \
                perl \
                wget 

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
        && locale-gen en_US.utf8 \
        && /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND noninteractive


COPY conda_dependencies.yml .
RUN mamba env create -f conda_dependencies.yml

RUN echo "source activate rnaseq-nf" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq-nf/bin:$PATH

# clean conda environment
RUN conda clean --yes --all
