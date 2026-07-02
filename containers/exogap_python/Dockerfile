FROM python:3.9.23-slim-bookworm

LABEL maintainer="Dorine MERLAT"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Dépendances système utiles pour certaines libs Python
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    git \
    wget \
    curl \
    libsqlite3-dev \
    && rm -rf /var/lib/apt/lists/*

# Mise à jour pip
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Librairies Python non standard nécessaires
RUN pip install --no-cache-dir \
    biopython \
    ete3 \
    gffpandas \
    gffutils \
    matplotlib \
    numpy \
    pandas \
    requests \
    seaborn

WORKDIR /app

CMD ["python"]