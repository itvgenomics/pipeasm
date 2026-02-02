FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    zlib1g-dev \
    libhts-dev \
    libbz2-dev \
    r-base \
    ca-certificates \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# ---------- FASTK ----------
RUN wget https://github.com/thegenemyers/FASTK/archive/refs/tags/v1.2.tar.gz \
    && tar -xzf v1.2.tar.gz \
    && rm v1.2.tar.gz

WORKDIR /opt/FASTK-1.2
RUN make

# ---------- MERQURY.FK ----------
WORKDIR /opt
RUN wget https://github.com/thegenemyers/MERQURY.FK/archive/refs/tags/v1.2.tar.gz \
    && tar -xzf v1.2.tar.gz \
    && rm v1.2.tar.gz

WORKDIR /opt/MERQURY.FK-1.2
RUN make

# ---------- GENESCOPE.FK ----------
WORKDIR /opt
RUN wget https://github.com/thegenemyers/GENESCOPE.FK/archive/refs/tags/v1.0.0.tar.gz \
    && tar -xzf v1.0.0.tar.gz \
    && rm v1.0.0.tar.gz

WORKDIR /opt/GENESCOPE.FK-1.0.0

# R library location
RUN mkdir -p /opt/R_libs

# Install R dependencies
RUN Rscript -e 'install.packages(c("argparse", "minpack.lm", "ggplot2", "viridis"), lib="/opt/R_libs", repos="https://cloud.r-project.org")'

# Run project install script
RUN sed -i 's|local_lib_path = "~/R_libs/"|local_lib_path = "/opt/R_libs/"|' /opt/GENESCOPE.FK-1.0.0/install.R
RUN Rscript install.R

# Set R library path
ENV R_LIBS=/opt/R_libs

# Set PATH
ENV PATH="/opt/GENESCOPE.FK-1.0.0:/opt/FASTK-1.2:/opt/MERQURY.FK-1.2:${PATH}"

# ---------- Smudgeplot ----------
WORKDIR /opt
RUN wget https://github.com/KamilSJaron/smudgeplot/archive/refs/tags/v0.5.3.tar.gz \
    && tar -xzf v0.5.3.tar.gz \
    && rm v0.5.3.tar.gz

WORKDIR /opt/smudgeplot-0.5.3
RUN make \
    && python3 -m pip install . --break-system-packages

WORKDIR /

CMD ["/bin/bash"]
