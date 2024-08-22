FROM python:3.12
LABEL authors="chris.cheshire@crick.ac.uk"

# Update pip to latest version
RUN python -m pip install --upgrade pip

# Add thesource files to the image
COPY . /usr/src/genome_tools
WORKDIR /usr/src/genome_tools

# Update version
RUN pip install toml
RUN python update_version.py

# Install program
RUN pip install .

CMD ["bash"]
