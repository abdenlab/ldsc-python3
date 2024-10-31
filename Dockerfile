# --------------------- Stage 1: Build the wheel ---------------------
# Use the official Ubuntu 23.04 as the base image for the builder
FROM ubuntu:23.04 AS builder

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    POETRY_VIRTUALENVS_CREATE=false \
    POETRY_NO_INTERACTION=1

# Update and install system dependencies for building
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    python3.11 \
    python3.11-venv \
    python3.11-dev \
    && rm -rf /var/lib/apt/lists/*

# Ensure python3 and pip3 are the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 10

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 -

# Add Poetry to PATH
ENV PATH="/root/.local/bin:$PATH"

# Set the working directory inside the container
WORKDIR /app

# Copy pyproject.toml and poetry.lock if available
COPY pyproject.toml poetry.lock* /app/

# Install project dependencies without installing the package itself
RUN poetry install --only main

# Copy the rest of the project files
COPY . /app

# Install the package
RUN poetry install

# Set the entrypoint to ldsc
ENTRYPOINT ["poetry", "run"]