# SynTerra Container
# All bioinformatics dependencies in one container

FROM mambaorg/micromamba:1.5.8

LABEL maintainer="SynTerra Developers"
LABEL description="Container for SynTerra synteny-guided gene finder"

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Create environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy bin scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/synterra/bin/
ENV PATH="/opt/synterra/bin:$PATH"

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
