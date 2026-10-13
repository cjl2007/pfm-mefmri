ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE}

SHELL ["/bin/bash", "-lc"]

WORKDIR /opt/mefmri
COPY . /opt/mefmri

# Runtime defaults. External dependencies (FSL/Workbench/FreeSurfer/SimNIBS/MCR)
# should be provided by BASE_IMAGE or derived images.
ENV MEDIR=/opt/mefmri \
    HCPPIPEDIR=/opt/mefmri/HCPpipelines-master

ENTRYPOINT ["/opt/mefmri/bin/mefmri_pipeline.sh"]
