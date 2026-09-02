FROM docker:29.2.0-cli@sha256:f7da93bb3a6c57789eec1180058c44127f70548977cc18e602d52082506de826 AS docker_cli

FROM mambaorg/micromamba:2.3.3@sha256:e583f8e151fd1ebbc78b05cf4feb37ddbb518892dbabf6c0c03ea8ad613bc519

ARG BUILD_VERSION=1.0.0
ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="lowpass-variants-nextflow" \
      org.opencontainers.image.description="Low-pass SNV/indel workflow runtime for fresh and FFPE samples" \
      org.opencontainers.image.source="https://github.com/cfarkas/lowpass-variants-nextflow" \
      org.opencontainers.image.version="${BUILD_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && rm -f /tmp/environment.yml

USER root
COPY --from=docker_cli /usr/local/bin/docker /usr/local/bin/docker
RUN chmod 0755 /usr/local/bin/docker

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/lowpass-variants-nextflow
USER $MAMBA_USER

ENV PATH=/opt/conda/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    NXF_DISABLE_CHECK_LATEST=true
WORKDIR /opt/lowpass-variants-nextflow
CMD ["/bin/bash"]
