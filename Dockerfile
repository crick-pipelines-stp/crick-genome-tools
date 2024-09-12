FROM python:3.12
LABEL authors="chris.cheshire@crick.ac.uk"

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Update pip to latest version
RUN python -m pip install --upgrade pip

# Add thesource files to the image
COPY . /usr/src/crick_genome_tools
WORKDIR /usr/src/crick_genome_tools

# Update version
RUN pip install toml
RUN python update_version.py

# Install program
RUN pip install .

CMD ["bash"]
