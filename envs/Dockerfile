FROM mambaorg/micromamba:2.5.0

# This PLINK build targets Terra's Linux x86-64 workers.
RUN micromamba install --yes --name base --strict-channel-priority \
      --channel conda-forge --channel bioconda \
      'python=3.11.11' \
      'plink2=2.0.0a.6.9=h9948957_0' \
      'zstandard=0.23.0' \
    && micromamba clean --all --yes

# Cromwell can override the image entrypoint. Keep tools on PATH.
ENV PATH="/opt/conda/bin:${PATH}"
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/pair_ld.py /opt/pair_ld/pair_ld.py
