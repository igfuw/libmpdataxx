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
        software-properties-common \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update -qq \
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
        gcc-6 g++-6 \
        python3-dev \
        python3-h5py \
        python3-scipy \
        python3-matplotlib \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ARG PYVER

RUN if [ "$PYVER" = 3 ]; then sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 10; fi

RUN wget -O /usr/local/include/gnuplot-iostream.h https://raw.githubusercontent.com/dstahlke/gnuplot-iostream/master/gnuplot-iostream.h

COPY . .

WORKDIR /usr/local/src/libmpdataxx/libmpdata++/build

RUN cmake .. \
    && make \
    && make install

WORKDIR /usr/local/src/libmpdataxx/
