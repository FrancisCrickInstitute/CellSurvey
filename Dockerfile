FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH="/root/.pixi/bin:${PATH}"

WORKDIR /app

# Pixi environment to build: "default" (Stardist/TensorFlow) or "instanseg" (PyTorch).
#   docker build -t cellsurvey .
#   docker build --build-arg PIXI_ENV=instanseg -t cellsurvey-instanseg .
ARG PIXI_ENV=default
ENV PIXI_ENV=${PIXI_ENV}

# Copy pixi config for dependency resolution
COPY pixi.toml ./

RUN pixi install -e ${PIXI_ENV}

# Copy the pipeline code
COPY run.py ./
COPY cellsurvey ./cellsurvey

ENV TF_USE_LEGACY_KERAS=1

ENTRYPOINT ["sh", "-c", "exec pixi run -e \"$PIXI_ENV\" python run.py \"$@\"", "--"]
