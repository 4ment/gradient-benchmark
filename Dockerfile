FROM continuumio/anaconda3
LABEL "author"="Mathieu Fourment"
LABEL "company"="University of Technology Sydney"

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
		autoconf \
		automake \
		build-essential \
		cmake \
		default-jdk \
		git \
		libgsl0-dev \
		libtool \
		pkg-config \
		python2.7 \
		python-dev \
		python-tk \
		unzip \
		wget \
		zlib1g-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/phylovi/bito /bito
WORKDIR /bito
RUN git submodule update --init --recursive
RUN conda env create -f environment.yml
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate bito && cd /bito && make
WORKDIR /

# Programs for treetime_validation
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py && python2.7 get-pip.py
RUN pip2 install phylo-treetime==0.7.4 click biopython==1.76

RUN wget http://www.microbesonline.org/fasttree/FastTree-2.1.11.c
RUN gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree-2.1.11.c -lm && mv FastTree /usr/local/bin

RUN wget https://github.com/tothuhien/lsd/archive/refs/heads/master.zip && unzip master.zip && chmod +x lsd-master/bin/lsd_unix && mv lsd-master/bin/lsd_unix /usr/local/bin/lsd

RUN wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.8.4/BEASTv1.8.4.tgz
RUN tar -xzvf BEASTv1.8.4.tgz && ln -s /BEASTv1.8.4/bin/beast /usr/local/bin

RUN wget https://github.com/tothuhien/lsd2/releases/download/v.2.3/lsd2_unix && chmod +x lsd2_unix && mv lsd2_unix /usr/local/bin/lsd2

# autodiff programs
RUN git clone --depth 1 https://github.com/4ment/physher.git /physher
WORKDIR /physher/Release
RUN cd /physher/Release && cmake -DBUILD_SHARED_LIBS=OFF -DBUILD_BENCHMARKING=1 .. && make && make install
RUN ln -s /physher/Release/examples/benchmarking /usr/local/bin/physher-benchmark
WORKDIR /

RUN git clone --depth 1 https://github.com/4ment/torchtree /torchtree
RUN cd /torchtree && /opt/conda/envs/bito/bin/pip install . && /opt/conda/envs/bito/bin/torchtree --help
RUN ln -s /torchtree/benchmarks/benchmark.py /usr/local/bin/torchtree-benchmark \
    && chmod +x /usr/local/bin/torchtree-benchmark

RUN git clone --depth 1 https://github.com/4ment/torchtree-bito /bitorch
RUN cd /bitorch && /opt/conda/envs/bito/bin/pip install .
RUN ln -s /bitorch/benchmarks/benchmark.py /usr/local/bin/bitorch-benchmark \
    && chmod +x /usr/local/bin/bitorch-benchmark
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate bito && /usr/local/bin/bitorch-benchmark --help

RUN git clone --depth 1 https://github.com/4ment/phylojax /phylojax
RUN cd /phylojax && /opt/conda/envs/bito/bin/pip install jax==0.2.24 jaxlib . && /opt/conda/envs/bito/bin/phylojax --help
RUN ln -s /phylojax/benchmarks/benchmark.py /usr/local/bin/phylojax-benchmark \
    && chmod +x /usr/local/bin/phylojax-benchmark

RUN git clone --depth 1 https://github.com/4ment/phylostan /phylostan
RUN cd /phylostan && pip install . && phylostan --help

RUN git clone --depth 1 https://github.com/christiaanjs/treeflow.git /treeflow
RUN cd /treeflow && /opt/conda/envs/bito/bin/pip install . 
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate bito && treeflow_benchmark --help

RUN echo "source activate bito" > ~/.bashrc
ENV PATH /opt/conda/envs/bito/bin:$PATH