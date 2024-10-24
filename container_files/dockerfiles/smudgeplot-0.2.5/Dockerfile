FROM continuumio/miniconda3:23.10.0-1
RUN conda install -c conda-forge -c bioconda smudgeplot=0.2.5 python=3.10 -y
WORKDIR /opt
RUN mkdir kmc
WORKDIR /opt/kmc
RUN wget https://github.com/refresh-bio/KMC/releases/download/v3.2.4/KMC3.2.4.linux.x64.tar.gz && tar xvzf *.tar.gz && rm *.tar.gz
ENV PATH="${PATH}:/opt/kmc/bin"
WORKDIR /