FROM nfcore/base:1.7
LABEL authors="Darcy Jones" \
      description="Docker image containing all requirements for nf-core/catastroflow pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-catastroflow-1.0dev/bin:$PATH
