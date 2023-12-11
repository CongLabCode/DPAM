ARG CUDA=11.1.1
FROM nvidia/cuda:11.1.1-cudnn8-runtime-ubuntu18.04
ARG CUDA


SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        build-essential \
        cmake \
        cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
        git \
        tzdata \
        wget \
        dialog \
        gfortran \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

RUN wget -q -P /tmp \
  https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh \
    && bash /tmp/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"
RUN conda install -y -c bioconda blast-legacy
RUN conda install -y -c biocore psipred
RUN conda install -y -c salilab dssp


RUN pip install numpy
RUN pip install tensorflow==1.14
RUN pip install protobuf==3.20.*


COPY utilities/DaliLite.v5.tar.gz /opt
RUN cd /opt \
    && tar -zxvf DaliLite.v5.tar.gz \
    && cd /opt/DaliLite.v5/bin \
    && make clean \
    && make \
    && ln -s /opt/DaliLite.v5/bin/dali.pl /usr/bin \
    && rm  /opt/DaliLite.v5.tar.gz

RUN git clone https://github.com/soedinglab/pdbx.git /opt/pdbx \
    && mkdir /opt/pdbx/build \
    && pushd /opt/pdbx/build \
    && cmake ../ \
    && make install \
    && popd

RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
    && mkdir /tmp/hh-suite/build \
    && pushd /tmp/hh-suite/build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
    && make -j 4 && make install \
    && ln -s /opt/hhsuite/bin/* /usr/bin \
    && popd \
    && rm -rf /tmp/hh-suite

RUN mkdir /opt/DPAM && mkdir /opt/DPAM/scripts
COPY scripts/*.py /opt/DPAM/scripts
COPY utilities/HHPaths.pm /opt/hhsuite/scripts
COPY utilities/pdb2fasta /usr/bin
COPY utilities/foldseek /usr/bin

RUN chmod -R +x /opt/DPAM/scripts

ENV PATH="/opt/DPAM/scripts:/opt/hhsuite/scripts:/opt/hhsuite/bin:/opt/DaliLite.v5/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"
ENV PERL5LIB="/usr/local/lib/perl5:/opt/hhsuite/scripts"
ENV OPENBLAS_NUM_THREADS=1

