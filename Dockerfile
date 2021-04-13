# Dockerfile for Simons Brain Aging Pheweb site
FROM ubuntu:18.04
MAINTAINER Kevin M Anderson "kevinanderson@fas.harvard.edu"

WORKDIR /app

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        sudo \
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

RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo
USER docker