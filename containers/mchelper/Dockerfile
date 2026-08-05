FROM continuumio/miniconda3

LABEL maintainer="Dorine MERLAT <dorine.merlat@etu.unistra.fr>"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /opt

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        libgl1 \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 && \
    rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/GonzalezLab/MCHelper.git /opt/MCHelper && \
    rm -rf /opt/MCHelper/.git

RUN conda env create -f /opt/MCHelper/MCHelper.yml && \
    conda clean -afy && \
    find /opt/conda -type d -name '__pycache__' -prune -exec rm -rf {} + && \
    find /opt/conda -type f \( -name '*.pyc' -o -name '*.pyo' -o -name '*.a' \) -delete && \
    rm -rf \
        /opt/conda/envs/MCHelper/include \
        /opt/conda/envs/MCHelper/share/doc \
        /opt/conda/envs/MCHelper/share/man \
        /opt/conda/envs/MCHelper/share/info

RUN conda run -n MCHelper python -c "import cv2; print('cv2 OK')"

RUN printf '#!/bin/sh\nexec conda run --no-capture-output -n MCHelper python /opt/MCHelper/MCHelper.py "$@"\n' \
    > /usr/local/bin/MCHelper.py && \
    chmod +x /usr/local/bin/MCHelper.py

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "MCHelper", "python", "/opt/MCHelper/MCHelper.py"]