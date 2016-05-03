FROM ubuntu:latest
MAINTAINER "zack ramjan" <zack.hr@gREMOVEmail.com>
ENV container docker
RUN apt-get update -y && \
    apt-get install git -y && \
    apt-get install build-essential -y && \
    apt-get install zlib1g-dev -y && \
    apt-get install libncurses5-dev -y && \
    apt-get install wget -y

RUN cd / && \
    git clone https://github.com/zwdzwd/biscuit.git && \
    cd /biscuit && \
    pwd && \
    make && \
    cp /biscuit/bin/biscuit /usr/bin && \
    cd / && \
    wget "https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2" && \
    tar -jxvf samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    make && \
    cp samtools /usr/bin 

VOLUME [ "/data" ]

#how i run:
#docker run -ti -d --name=biscuit_docker -v /root/biscuit:/data biscuit

#connect to live instance
#docker exec -it biscuit_docker /bin/bash

