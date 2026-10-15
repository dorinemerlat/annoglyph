FROM python:3.9.23-slim-bookworm

LABEL maintainer="Dorine MERLAT"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install system dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    git \
    wget \
    curl \
    procps \
    libsqlite3-dev \
    && rm -rf /var/lib/apt/lists/*

# Upgrade Python packaging tools.
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel

COPY requirements.txt .

RUN python -m pip install --no-cache-dir -r requirements.txt

# Verify Python dependencies
RUN python -c "import Bio, ete3, gffpandas, gffutils, matplotlib, numpy, pandas, requests, seaborn"

WORKDIR /app

CMD ["python"]