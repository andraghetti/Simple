# Simple Pipeline Docker Container - Multi-Architecture Support
# Supports both AMD64 and ARM64 architectures

FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV MAKEFLAGS="-j$(nproc)"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    openjdk-8-jdk \
    wget \
    curl \
    build-essential \
    software-properties-common \
    apt-transport-https \
    ca-certificates \
    gnupg \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R and packages
RUN apt-get update && apt-get install -y \
    r-base \
    r-cran-ggplot2 \
    r-cran-reshape2 \
    && rm -rf /var/lib/apt/lists/*

# Set Java environment (auto-detect architecture)
RUN ARCH=$(dpkg --print-architecture) && \
    echo "Detected architecture: $ARCH" && \
    if [ "$ARCH" = "arm64" ]; then \
        echo "export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-arm64" >> /etc/environment; \
    else \
        echo "export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64" >> /etc/environment; \
    fi && \
    echo "export PATH=\$JAVA_HOME/bin:\$PATH" >> /etc/environment

# Set default Java environment
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# Verify R packages are available
RUN R --slave -e "library(ggplot2); library(reshape2); cat('R packages loaded successfully\n')"

# Set working directory
WORKDIR /app

# Copy all files (excluding output and data directories via .dockerignore)
COPY . /app/

# Create runs directory for output
RUN mkdir -p /app/runs /app/refs

# Install BWA 0.7.19 (latest bug-fixed version) with architecture-specific optimizations
RUN echo "Installing BWA 0.7.19..." && \
    wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.19.tar.gz && \
    tar -xzf v0.7.19.tar.gz && \
    cd bwa-0.7.19 && \
    ARCH=$(dpkg --print-architecture) && \
    if [ "$ARCH" = "arm64" ]; then \
        make CFLAGS="-O2 -march=armv8-a"; \
    else \
        make CFLAGS="-O2 -march=x86-64"; \
    fi && \
    cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa-0.7.19*

# Install latest samtools (1.22.1)
RUN echo "Installing samtools 1.22.1..." && \
    wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2 && \
    tar -xjf samtools-1.22.1.tar.bz2 && \
    cd samtools-1.22.1 && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf samtools-1.22.1*

# Install latest Picard (2.27.5)
RUN echo "Installing Picard 2.27.5..." && \
    wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -O /usr/local/bin/picard.jar

# Keep existing GATK 3.7.0 for compatibility (already in programs folder)
# We'll copy it to system location for consistency
RUN echo "Setting up GATK..." && \
    cp programs/GenomeAnalysisTK.jar /usr/local/bin/GenomeAnalysisTK.jar

# Keep existing snpEff 4.3p for compatibility (already in programs folder)
# We'll copy it to system location for consistency
RUN echo "Setting up snpEff..." && \
    cp -r programs/snpEff /usr/local/ && \
    ln -s /usr/local/snpEff/snpEff.jar /usr/local/bin/snpEff.jar && \
    ln -s /usr/local/snpEff/SnpSift.jar /usr/local/bin/SnpSift.jar

# Make scripts executable
RUN chmod +x ./scripts/simple.sh

# Set entrypoint
ENTRYPOINT ["./scripts/simple.sh"]
