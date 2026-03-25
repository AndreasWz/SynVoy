# SynVoy Container
# All bioinformatics dependencies in one container

FROM mambaorg/micromamba:1.5.8

LABEL maintainer="SynVoy Developers"
LABEL description="Container for SynVoy synteny-guided gene finder"

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Create environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy bin scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/synvoy/bin/
ENV PATH="/opt/synvoy/bin:$PATH"

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
