################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################
FROM continuumio/miniconda3
MAINTAINER Caspar Gross <mail@caspar.one>
LABEL description="contains all the dependencies for plasmid Identification pipeline at github.com/caspargross/plasmidIdentification" 

SHELL ["/bin/bash", "-c"]

# Install conda envs
ADD env/PI_env.yml /tmp/PI_env.yml
RUN conda env create -f /tmp/PI_env.yml -q && conda clean -a

# Download CARD-Antibiotic resistance database
RUN wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2 && bash && source activate PI_env && rgi load --afile card.json
