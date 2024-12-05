# Base stage: Common setup for both test and prod
FROM condaforge/miniforge3

# Install apt packages
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps \
    locales \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Generate and set the locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Set up conda channels with priority
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# Install conda packages
RUN mamba install -y python=3.12
RUN mamba install -y bwa-mem2
RUN mamba install -y samtools
RUN mamba install -y bcftools
RUN mamba install -y abra2
RUN mamba install -y pysam
RUN mamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /home

# Copy the source files
COPY . /usr/src/crick_genome_tools
WORKDIR /usr/src/crick_genome_tools

# Install the program
RUN pip install pytest
RUN pip install .

CMD ["bash"]
