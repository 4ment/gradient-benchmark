FROM continuumio/anaconda3:2022.10
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

RUN git clone https://github.com/phylovi/bito /bito
WORKDIR /bito
RUN git checkout cc0806abcd0b9f2fab604e800c674c9a5c5afebe
RUN git submodule update --init --recursive
RUN sed -i 's/>= 3.7/== 3.10/g' environment.yml
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

RUN git clone https://github.com/4ment/physher.git /physher
RUN cd /physher && git checkout b19ff2f9422f29ba1ab31306a3fe29ab6a6f607b
WORKDIR /physher/Release
RUN cd /physher/Release && cmake -DBUILD_SHARED_LIBS=OFF -DBUILD_BENCHMARKING=1 .. && make && make install
RUN ln -s /physher/Release/examples/benchmarking /usr/local/bin/physher-benchmark
WORKDIR /

# autodiff programs
RUN git clone https://github.com/4ment/torchtree /torchtree
RUN cd /torchtree && git checkout f3831650a807e74cc2e9478009e57a41f47bed8d \
    && /opt/conda/envs/bito/bin/pip install torch==1.12.1 numpy==1.22 . \
	&& /opt/conda/envs/bito/bin/torchtree --help
RUN ln -s /torchtree/benchmarks/benchmark.py /usr/local/bin/torchtree-benchmark \
    && chmod +x /usr/local/bin/torchtree-benchmark

RUN git clone https://github.com/4ment/torchtree-bito /bitorch
RUN cd /bitorch && git checkout e2a95cefb13968f95f6e5520bd0a52d726ee7fc9 \
    && /opt/conda/envs/bito/bin/pip install .
RUN ln -s /bitorch/benchmarks/benchmark.py /usr/local/bin/bitorch-benchmark \
    && chmod +x /usr/local/bin/bitorch-benchmark
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate bito && /usr/local/bin/bitorch-benchmark --help

RUN git clone https://github.com/4ment/phylojax /phylojax
RUN cd /phylojax && git checkout a1612cae36292af76e8d24cc40d6544162c987aa \
    && /opt/conda/envs/bito/bin/pip install jax==0.2.24 jaxlib==0.3.7 numpy==1.22 . \
	&& /opt/conda/envs/bito/bin/phylojax --help
RUN ln -s /phylojax/benchmarks/benchmark.py /usr/local/bin/phylojax-benchmark \
    && chmod +x /usr/local/bin/phylojax-benchmark

RUN pip install phylostan==1.0.5 pystan==2.19.1.1 && phylostan --help

RUN git clone https://github.com/christiaanjs/treeflow.git /treeflow
RUN cd /treeflow && git checkout e3414dcc9e764d06abc3e19c1d0f55110499e2ea \
    && /opt/conda/envs/bito/bin/pip install tensorflow==2.10.0 tensorflow-probability==0.18.0 . 
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate bito && treeflow_benchmark --help

RUN echo "source activate bito" > ~/.bashrc
ENV PATH /opt/conda/envs/bito/bin:$PATH
