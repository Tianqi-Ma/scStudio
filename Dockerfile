# scStudio — dependency-free distribution.
# Everything (R + Seurat + Bioconductor + the app) is baked into this image, so a
# user only needs Docker. Build once, then anyone can run:
#
#   docker run --rm -p 3838:3838 -m 16g <image>
#   # then open http://localhost:3838
#
# The Bioconductor base image ships R plus a prebuilt Bioconductor toolchain,
# which makes the Bioc dependencies (scDblFinder, SingleR, scater/scran) install
# quickly and reliably.

FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL org.opencontainers.image.title="scStudio" \
      org.opencontainers.image.description="Local interactive single-cell RNA-seq analysis app" \
      org.opencontainers.image.source="https://github.com/Tianqi-Ma/scStudio"

# System libs occasionally needed by leiden/igraph/plotly stacks are already in
# the Bioconductor base. Install R package dependencies in a cached layer.
RUN R -e "install.packages(c( \
      'shiny','bslib','ggplot2','Matrix','plotly','DT','shinyWidgets', \
      'promises','future','progressr','remotes','Seurat','SeuratObject', \
      'harmony','clustree','ggrastr','scattermore'), \
      repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
      'SingleCellExperiment','SummarizedExperiment','scater','scran', \
      'scDblFinder','glmGamPoi','SingleR','celldex'), update=FALSE, ask=FALSE)"

# Install the app itself (copy source and install from local path).
WORKDIR /opt/scStudio
COPY . /opt/scStudio
RUN R -e "remotes::install_local('/opt/scStudio', dependencies = FALSE, upgrade = 'never')"

EXPOSE 3838

# Bind to 0.0.0.0 so the host browser can reach the container; do NOT auto-open a
# browser inside the container. Users open http://localhost:3838 themselves.
CMD ["R", "-e", "scStudio::run_app(host='0.0.0.0', port=3838, launch.browser=FALSE)"]
