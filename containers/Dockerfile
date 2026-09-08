# Docker image for Vedette. Build from the repository root:
#   docker build -f containers/Dockerfile -t vedette:latest .
FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md LICENSE /opt/vedette/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/vedette/src
RUN micromamba run -n base pip install --no-cache-dir /opt/vedette

ENTRYPOINT ["micromamba", "run", "-n", "base", "vedette"]
CMD ["--help"]
