FROM ubuntu:16.04

ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /usr/local/src/libmpdataxx

RUN apt-get update -qq \
    && apt-get install -yq --no-install-recommends \
        sudo \
        apt-utils \
        build-essential \
        pkg-config \
        git \
        cmake \
        ca-certificates \
        clang-4.0 \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update -qq \
    && apt-get install -yq --no-install-recommends \
        libblitz0-dev \
        libboost-all-dev \
        gnuplot-nox \
        wget \
        libhdf5-dev \
        hdf5-tools \
        python-h5py \
        python-scipy \
        python-matplotlib \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget -O /usr/local/include/gnuplot-iostream.h https://raw.githubusercontent.com/dstahlke/gnuplot-iostream/master/gnuplot-iostream.h

COPY . .

WORKDIR /usr/local/src/libmpdataxx/libmpdata++/build

RUN cmake .. && make

WORKDIR /usr/local/src/libmpdataxx/
