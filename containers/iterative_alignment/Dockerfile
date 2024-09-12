FROM mambaorg/micromamba:latest

# Install apt packages
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install conda packages
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER containers/iterative_alignment/env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV PATH /opt/conda/bin:$PATH
WORKDIR /home/mambauser
USER root

# Update pip to latest version
RUN python -m pip install --upgrade pip

# Add thesource files to the image
COPY . /usr/src/crick_genome_tools
WORKDIR /usr/src/crick_genome_tools

# Install program
RUN pip install '.[tests]'
USER $MAMBA_USER

CMD ["bash"]
