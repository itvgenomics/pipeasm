FROM continuumio/miniconda3:23.10.0-1
RUN conda install -c conda-forge -c bioconda merqury=1.3 openjdk=11 -y
RUN sed -i 's/`link //g' /opt/conda/bin/merqury.sh
RUN sed -i 's/`//g' /opt/conda/bin/merqury.sh