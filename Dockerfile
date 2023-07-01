FROM mambaorg/micromamba:1.4.4
MAINTAINER Reza Shahneh "mzare008@ucr.edu"

COPY conda-env.yml .
RUN micromamba env create -f conda-env.yml --name mod-site
RUN bash -c "source activate mod-site && pip install gunicorn"

COPY . /app
WORKDIR /app
