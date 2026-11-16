# Docker image for SVengers. Build from the repository root:
#   docker build -f containers/Dockerfile -t svengers:latest .
FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md LICENSE /opt/svengers/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/svengers/src
RUN micromamba run -n base pip install --no-cache-dir /opt/svengers

ENTRYPOINT ["micromamba", "run", "-n", "base", "svengers"]
CMD ["--help"]
