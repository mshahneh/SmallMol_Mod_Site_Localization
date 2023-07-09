<<<<<<< HEAD
FROM continuumio/miniconda3:4.10.3
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

RUN apt-get update && apt-get install -y build-essential
RUN conda install -c conda-forge mamba
COPY conda-env.yml .
RUN mamba env create -f conda-env.yml -n mod-site

RUN /bin/bash -c 'source activate mod-site'
RUN /bin/bash -c 'mamba install -c conda-forge gunicorn'
=======
FROM mambaorg/micromamba:1.4.4
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

COPY conda-env.yml .
RUN micromamba env create -f conda-env.yml --name mod-site
RUN echo "source activate mod-site" > ~/.bashrc
ENV PATH /opt/conda/envs/mod-site/bin:$PATH
>>>>>>> 6172b69 (chaning the package manager to mamba)

COPY . /app
WORKDIR /app
