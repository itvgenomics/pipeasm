FROM ubuntu:22.04
RUN apt update && apt install -y wget gcc xz-utils
WORKDIR /opt
RUN wget https://github.com/marbl/meryl/releases/download/v1.3/meryl-1.3.Linux-amd64.tar.xz
RUN tar -xJf meryl-1.3.*.tar.xz && rm *.tar.xz
ENV PATH="${PATH}:/opt/meryl-1.3/bin"
