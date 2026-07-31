FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH="/root/.pixi/bin:${PATH}"

WORKDIR /app

COPY pixi.toml pixi.lock latest.zip ./

RUN pixi install

ENV POT_BACKEND=numpy
ENV TF_USE_LEGACY_KERAS=1

ENTRYPOINT ["pixi", "run", "python", "run.py"]
