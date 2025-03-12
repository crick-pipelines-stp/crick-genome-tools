FROM condaforge/miniforge3

# Set up conda channels with priority
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# Install conda packages
RUN mamba install -y python=3.12
RUN mamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /home

# Copy the source files
COPY . crick_genome_tools
WORKDIR crick_genome_tools

# Install the program
RUN pip install .

# Expose Streamlit's default port
EXPOSE 8501

CMD ["bash"]
