FROM ubuntu:22.04
RUN apt update && apt install -y wget build-essential python3 python-is-python3 libz-dev
WORKDIR /opt
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.20.0.tar.gz
RUN tar -xvzf 0.20.0.tar.gz
RUN rm 0.20.0.tar.gz
WORKDIR /opt/hifiasm-0.20.0
RUN make
ENV PATH="${PATH}:/opt/hifiasm-0.20.0"

