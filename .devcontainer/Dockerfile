# Note: You can use any Debian/Ubuntu based image you want. 
FROM mcr.microsoft.com/devcontainers/base:bullseye

WORKDIR /tmp

RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.4.tar.gz
RUN tar -zxvf ViennaRNA-2.6.4.tar.gz
RUN cd ViennaRNA-2.6.4 && ./configure && make && sudo make install

# [Optional] Uncomment this section to install additional OS packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends ncbi-blast+