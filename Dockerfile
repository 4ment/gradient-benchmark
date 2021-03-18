FROM quay.io/matsengrp/conda-beagle

RUN apt-get update && apt-get install -y --no-install-recommends \
	cmake \
	libgsl0-dev \
	unzip \
	wget \
	git

RUN git clone --depth 1 --branch=main https://github.com/phylovi/libsbn /libsbn
WORKDIR /libsbn
RUN git submodule update --init --recursive

RUN /opt/conda/bin/conda env create -f environment.yml
RUN echo "conda activate libsbn" >> $HOME/.bashrc
RUN /opt/conda/envs/libsbn/bin/pip install phylostan

WORKDIR /
ENV BEAGLE_PREFIX /usr/local
ENV LD_LIBRARY_PATH /usr/local/lib

RUN wget https://github.com/4ment/physher/archive/marginal-v1.1.zip && unzip marginal-v1.1.zip
WORKDIR /physher-marginal-v1.1/Release
RUN cmake -DBUILD_SHARED_LIBS=OFF .. && make && make install

WORKDIR /data