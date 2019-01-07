
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# transmission: A tool for inferring endosymbiont biology from metagenome data.
# Copyright (C) 4-11-2018 Mark Juers

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

ARG BASE_CONTAINER=jupyter/r-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Mark Juers <mpjuers@indiana.edu>"

RUN mkdir Transmission
COPY --chown=1000 ./ Transmission/

USER root

ENV DEBIAN_FRONTEND noninteractive
ENV PATH "/home/jovyan/.local/bin:${PATH}"
RUN apt-get update && apt-get -yq install libreadline-dev

USER $NB_UID

RUN conda install pip
# Install build requirements.
RUN cd Transmission && \
    while read requirement; do \
        if ! conda install --yes $requirement; then \
            pip install --user $requirement; fi; done < requirements.txt
# Install R dependencies.
RUN Rscript -e 'install.packages(c("abc", "kdensity"), \
        repos = "https://ftp.ussg.iu.edu/CRAN/")'

RUN cd Transmission && python setup.py sdist && pip install --user .

