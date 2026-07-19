# metaQII-nf runtime image
# -----------------------------------------------------------------------------
# Bakes the pinned QIIME 2 amplicon environment (plus fastqc / itsxpress /
# q2-itsxpress) into a single image so runs are reproducible regardless of how
# conda channels drift over time.
#
# Build locally with Docker:
#   docker build -t ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0 .
#
# Push (optional) so HPC nodes can pull it:
#   docker push ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0
#
# Pull with Apptainer/Singularity on an HPC login node:
#   apptainer pull metaqii-nf_2026.1.0.sif \
#     docker://ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0
#
# Air-gapped alternative (no registry):
#   docker save ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0 -o metaqii.tar
#   # copy metaqii.tar to the cluster, then:
#   apptainer build metaqii-nf_2026.1.0.sif docker-archive://metaqii.tar
# -----------------------------------------------------------------------------
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="metaQII-nf" \
      org.opencontainers.image.description="QIIME 2 amplicon pipeline environment (16S & ITS)" \
      org.opencontainers.image.source="https://github.com/dariusng28-collab/metaQII-nf" \
      org.opencontainers.image.licenses="MIT"

# `ps` (procps) is required by Nextflow to collect per-task resource metrics.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Solve the pinned environment into the base env, then strip caches to slim the
# image. `q2-itsxpress` is installed explicitly with pip because `micromamba
# install -f` does not process the env file's `pip:` subsection.
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/qiime2-amplicon-2026.1-metaqii.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
    && micromamba run -n base pip install --no-cache-dir q2-itsxpress multiqc \
    && micromamba clean --all --yes \
    && rm -f /tmp/env.yml

# Put the base environment on PATH so `qiime`, `fastqc`, etc. are callable
# without an activation step (Nextflow invokes the process script directly).
ENV PATH="/opt/conda/bin:$PATH" \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

SHELL ["/bin/bash", "-c"]
