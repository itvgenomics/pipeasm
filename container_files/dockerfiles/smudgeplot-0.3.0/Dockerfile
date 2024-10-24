FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install --no-install-recommends -y \
    git build-essential libz-dev r-base libcurl4-gnutls-dev libxml2-dev libssl-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libbz2-dev libclang-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev
WORKDIR /opt
RUN git clone -b sploidyplot https://github.com/KamilSJaron/smudgeplot
WORKDIR /opt/smudgeplot
RUN Rscript -e 'install.packages(c("viridis", "argparse"))'
RUN Rscript -e 'install.packages("devtools")'
RUN make
WORKDIR /opt
RUN git clone https://github.com/thegenemyers/FASTK.git
WORKDIR /opt/FASTK
RUN make
ENV PATH=$PATH:"/opt/FASTK"
WORKDIR /