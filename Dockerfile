FROM nfcore/base:1.7
LABEL authors="Darcy Jones" \
      description="Docker image containing all requirements for catastrophy-pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a --yes
ENV PATH /opt/conda/envs/catastrophy-pipeline/bin:$PATH
