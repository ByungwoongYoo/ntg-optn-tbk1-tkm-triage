FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER lantern_submission/configs/environment-lock.yml /tmp/environment-lock.yml
RUN micromamba install --yes --name base --file /tmp/environment-lock.yml \
    && micromamba clean --all --yes

WORKDIR /opt/lantern
COPY --chown=$MAMBA_USER:$MAMBA_USER lantern_submission /opt/lantern/lantern_submission
COPY --chown=$MAMBA_USER:$MAMBA_USER lantern_v6 /opt/lantern/lantern_v6
COPY --chown=$MAMBA_USER:$MAMBA_USER lantern_cami3/scripts /opt/lantern/lantern_cami3/scripts
RUN chmod +x /opt/lantern/lantern_submission/run_frozen_v7.sh \
    /opt/lantern/lantern_submission/scripts/*.py

ENV PYTHONPATH=/opt/lantern/lantern_cami3/scripts:/opt/lantern/lantern_v6:/opt/lantern/lantern_submission/scripts
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["python", "lantern_submission/scripts/load_explicit_mapping.py", "--help"]
