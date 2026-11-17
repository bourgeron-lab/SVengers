# Docker image for SVengers. Build from the repository root:
#   docker build -f containers/Dockerfile -t svengers:latest .
# Published for every release as ghcr.io/bourgeron-lab/svengers:<version>.
# The image carries the filter model (/opt/svengers/models/v12b) and the generic
# GRCh38 config (/opt/svengers/configs/grch38.toml); the reference files are mounted
# from the directory built by workflows/build_grch38_resources.sh (docs/RESOURCES.md).
FROM mambaorg/micromamba:1.5.8
LABEL org.opencontainers.image.source="https://github.com/bourgeron-lab/SVengers"
LABEL org.opencontainers.image.description="SVengers: multi-evidence germline CNV caller (DEL/DUP)"
LABEL org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md LICENSE /opt/svengers/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/svengers/src
COPY --chown=$MAMBA_USER:$MAMBA_USER models/v12b /opt/svengers/models/v12b
COPY --chown=$MAMBA_USER:$MAMBA_USER configs/grch38.toml /opt/svengers/configs/grch38.toml
RUN micromamba run -n base pip install --no-cache-dir /opt/svengers

ENTRYPOINT ["micromamba", "run", "-n", "base", "svengers"]
CMD ["--help"]
