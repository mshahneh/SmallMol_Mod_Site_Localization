<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 0b9606e (dash app container)
FROM continuumio/miniconda3:4.10.3
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

RUN apt-get update && apt-get install -y build-essential
RUN conda install -c conda-forge mamba
COPY conda-env.yml .
RUN mamba env create -f conda-env.yml -n mod-site

RUN /bin/bash -c 'source activate mod-site'
RUN /bin/bash -c 'mamba install -c conda-forge gunicorn'
<<<<<<< HEAD
=======
FROM mambaorg/micromamba:1.4.4
=======
FROM continuumio/miniconda3:4.10.3
>>>>>>> 0b9606e (dash app container)
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

RUN apt-get update && apt-get install -y build-essential
RUN conda install -c conda-forge mamba
COPY conda-env.yml .
<<<<<<< HEAD
RUN micromamba env create -f conda-env.yml --name mod-site
<<<<<<< HEAD
RUN echo "source activate mod-site" > ~/.bashrc
ENV PATH /opt/conda/envs/mod-site/bin:$PATH
>>>>>>> 6172b69 (chaning the package manager to mamba)
=======
RUN bash -c "source activate mod-site && pip install gunicorn"
>>>>>>> 0cfd125 (added docker)
=======
RUN mamba env create -f conda-env.yml -n mod-site

RUN /bin/bash -c 'source activate mod-site'
RUN /bin/bash -c 'mamba install -c conda-forge gunicorn'
>>>>>>> 0b9606e (dash app container)
=======
>>>>>>> 0b9606e (dash app container)

COPY . /app
WORKDIR /app
