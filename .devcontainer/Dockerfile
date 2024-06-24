FROM mcr.microsoft.com/devcontainers/python:1-3.12-bullseye

# RUN apt update -y && apt install --no-install-recommends 

# Install ViennaRNA package dependency to running AddGenome.py

WORKDIR /tmp

RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.4.tar.gz
RUN tar -zxvf ViennaRNA-2.6.4.tar.gz
RUN cd ViennaRNA-2.6.4
RUN ./configure
RUN make
RUN sudo make install