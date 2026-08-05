FROM continuumio/miniconda3

LABEL maintainer="Dorine MERLAT <dorine.merlat@etu.unistra.fr>"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

ARG MCHELPER_REF=v1.7.0

WORKDIR /opt

# Install system dependencies.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        libgl1 \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        unzip && \
    rm -rf /var/lib/apt/lists/*

# Download MCHelper.
RUN git clone --depth 1 --branch "${MCHELPER_REF}" \
        https://github.com/GonzalezLab/MCHelper.git \
        /opt/MCHelper && \
    rm -rf /opt/MCHelper/.git

# Create the MCHelper Conda environment.
RUN conda env create -f /opt/MCHelper/MCHelper.yml && \
    conda clean -afy

ENV CONDA_PREFIX=/opt/conda/envs/MCHelper
ENV CONDA_DEFAULT_ENV=MCHelper
ENV PATH=/opt/conda/envs/MCHelper/bin:${PATH}

# Extract the databases distributed with MCHelper and create the BLAST database.
RUN cd /opt/MCHelper/db && \
    for archive in ./*.zip; do \
        [ -e "${archive}" ] || continue; \
        unzip -o "${archive}"; \
    done && \
    rm -f ./*.zip && \
    test -s allDatabases.clustered_rename.fa && \
    makeblastdb \
        -in allDatabases.clustered_rename.fa \
        -dbtype nucl

# Add the Pfam profiles used by MCHelper.
# The Pfam profile archive must be downloaded manually because the URGI website
# is protected by an anti-bot mechanism that prevents automated downloads.
#
# Download:
# https://urgi.versailles.inrae.fr/download/repet/profiles/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz
#
# Place the downloaded file next to this Dockerfile before building the image.
COPY ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz /tmp/

RUN cd /opt/MCHelper/db && \
    tar -tzf /tmp/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz >/dev/null && \
    tar -xzf /tmp/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz && \
    test -s ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm && \
    mv ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm Pfam35.0.hmm && \
    rm /tmp/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz

# Remove files that are not required at runtime.
RUN find /opt/conda -type d -name '__pycache__' -prune -exec rm -rf {} + && \
    find /opt/conda -type f \( -name '*.pyc' -o -name '*.pyo' -o -name '*.a' \) -delete && \
    rm -rf \
        /opt/conda/envs/MCHelper/include \
        /opt/conda/envs/MCHelper/share/doc \
        /opt/conda/envs/MCHelper/share/man \
        /opt/conda/envs/MCHelper/share/info

# Provide a convenient executable without forcing a container entrypoint.
RUN printf '%s\n' \
        '#!/bin/sh' \
        'exec python /opt/MCHelper/MCHelper.py "$@"' \
        > /usr/local/bin/MCHelper.py && \
    chmod +x /usr/local/bin/MCHelper.py

# Validate the main Python dependencies, executables and databases.
RUN python -c "import Bio, cv2, numpy, pandas, scipy, sklearn; print('Python dependencies OK')" && \
    command -v blastn && \
    command -v makeblastdb && \
    command -v hmmscan && \
    command -v hmmpress && \
    test -s /opt/MCHelper/db/allDatabases.clustered_rename.fa && \
    test -s /opt/MCHelper/db/allDatabases.clustered_rename.fa.nhr && \
    test -s /opt/MCHelper/db/allDatabases.clustered_rename.fa.nin && \
    test -s /opt/MCHelper/db/allDatabases.clustered_rename.fa.nsq && \
    test -s /opt/MCHelper/db/Pfam35.0.hmm && \
    MCHelper.py --help >/dev/null

WORKDIR /opt/MCHelper

CMD ["bash"]