FROM continuumio/miniconda3:4.10.3
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

RUN apt-get update && apt-get install -y build-essential
RUN conda install -c conda-forge mamba
COPY conda-env.yml .
RUN mamba env create -f conda-env.yml -n mod-site

RUN /bin/bash -c 'source activate mod-site'
RUN /bin/bash -c 'mamba install -c conda-forge gunicorn'

COPY . /app
WORKDIR /app
