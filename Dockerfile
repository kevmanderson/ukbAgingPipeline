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


RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'


ENV R_BASE_VERSION 3.6.3
ENV DEBIAN_FRONTEND=noninteractive

# Now install R and littler, and create a link for littler in /usr/local/bin
# Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		r-base=${R_BASE_VERSION}*
RUN apt-get install -y r-base-core=${R_BASE_VERSION}*
RUN apt-get install -y r-base-dev=${R_BASE_VERSION}*


# R
RUN R -e "install.packages('jsonlite',dependencies=TRUE, repos='http://cran.rstudio.com/')"



# miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py39_4.9.2-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py39_4.9.2-Linux-x86_64.sh



# Python modules
RUN python3 -m pip install --upgrade pip setuptools
RUN python3 -m pip install pandas==1.2.2
RUN python3 -m pip install numpy==1.20.1
RUN python3 -m pip install scipy==1.6.1
RUN python3 -m pip install numpy==1.20.1
RUN python3 -m pip install requests==2.25.1
RUN python3 -m pip install wget~=3.2
RUN python3 -m pip install pybiomart==0.2.0
RUN python3 -m pip install pysqlite3==0.4.6
RUN python3 -m pip install argparse==1.1



# plink 1.9
RUN mkdir /opt/plink \
&& cd /opt/plink \
&& wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip \
&& unzip plink_linux_x86_64_20201019.zip
ENV PATH /opt/plink:$PATH


# plink 2.0
RUN mkdir /opt/plink2 \
&& cd /opt/plink2 \
&& wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip \
&& unzip plink2_linux_avx2.zip
ENV PATH /opt/plink2:$PATH


# UK Biobank utilities
RUN wget -nd biobank.ndph.ox.ac.uk/showcase/util/ukbmd5
RUN wget -nd biobank.ndph.ox.ac.uk/showcase/util/ukbconv
RUN wget -nd biobank.ndph.ox.ac.uk/showcase/util/ukbunpack
RUN wget -nd biobank.ndph.ox.ac.uk/showcase/util/ukbfetch
RUN chmod 755 /app/*
ENV PATH /app:$PATH


RUN R -e "install.packages('jsonlite',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('feather',dependencies=TRUE, repos='http://cran.rstudio.com/')"


# copy all the scripts into root dir within the container
COPY scripts /scripts/
COPY ref_files /ref_files/






