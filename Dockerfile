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

# The tools' folder and the variables that activating the environment sets, as image variables:
# Docker runs the image's entrypoint, which activates the environment, but Singularity and
# Apptainer do not run it (singularity exec), and they keep these variables. Augustus reads its
# gene models from AUGUSTUS_CONFIG_PATH.
ENV PATH="/opt/conda/bin:/opt/synvoy/bin:$PATH" \
    CONDA_PREFIX="/opt/conda" \
    AUGUSTUS_CONFIG_PATH="/opt/conda/config/" \
    AUGUSTUS_SCRIPTS_PATH="/opt/conda/bin/" \
    AUGUSTUS_BIN_PATH="/opt/conda/bin/"

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
