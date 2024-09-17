# Base stage: Common setup for both test and prod
FROM mambaorg/micromamba:latest AS base

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

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /home/mambauser
USER root

# Update pip to the latest version
RUN python -m pip install --upgrade pip

# Production stage: Copies source code and installs it
FROM base AS prod

# Copy the source files
COPY . /usr/src/crick_genome_tools
WORKDIR /usr/src/crick_genome_tools

# Install the program
RUN pip install '.'

USER $MAMBA_USER
CMD ["bash"]

# Test stage: Installs pytest but does not copy source code
FROM base AS test

# Install pytest
RUN pip install pytest

USER $MAMBA_USER
CMD ["bash"]
