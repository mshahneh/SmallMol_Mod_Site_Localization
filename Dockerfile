FROM mambaorg/micromamba:1.4.4
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

COPY conda-env.yml .
RUN micromamba env create -f conda-env.yml --name mod-site
RUN echo "source activate mod-site" > ~/.bashrc
ENV PATH /opt/conda/envs/mod-site/bin:$PATH

COPY . /app
WORKDIR /app
