# Use the official Ubuntu 24.04 LTS as the base image
FROM ubuntu:24.04

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    POETRY_VIRTUALENVS_CREATE=false \
    POETRY_NO_INTERACTION=1

# Update and install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    wget \
    git \
    python3.12 \
    python3.12-venv \
    python3.12-dev \
    samtools \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Ensure python3 and pip3 are the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1 && \
    update-alternatives --install /usr/bin/pip3 pip3 /usr/bin/pip3.12 1

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 -

# Add Poetry to PATH
ENV PATH="/root/.local/bin:$PATH"

# Set the working directory inside the container
WORKDIR /app

# Copy pyproject.toml and poetry.lock if available
COPY pyproject.toml poetry.lock* /app/

# Install project dependencies
RUN poetry install --no-root --only main

# Copy the rest of the project files
COPY . /app

# Install the project
RUN poetry install --no-dev