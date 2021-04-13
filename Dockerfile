# Dockerfile for Simons Brain Aging Pheweb site
FROM ubuntu:18.04
MAINTAINER Kevin M Anderson "kevinanderson@fas.harvard.edu"

WORKDIR /app

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        libhdf5-serial-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# miniconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py39_4.9.2-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py39_4.9.2-Linux-x86_64.sh
ENV PATH /miniconda/bin:${PATH}
