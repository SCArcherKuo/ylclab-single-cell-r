FROM rocker/r-ver:4.4.0

# Layer 1: system libs (rarely changes)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# Layer 2: heavy R deps (rarely changes — cache this layer in CI)
# Use cloud.r-project.org for CRAN packages to get latest versions;
# matrixStats >= 1.4.1 is required by BiocManager's MatrixGenerics.
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "install.packages('matrixStats', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(c('scran', 'scater', 'glmGamPoi', 'MAST', 'DESeq2', 'ProteoMM'), ask=FALSE)" && \
    Rscript -e "install.packages(c('Seurat','harmony','optparse','jsonlite','arrow'), \
                                  repos='https://cloud.r-project.org')"

# Layer 3: ylclabscm package (rebuilt on each release)
COPY . /build/ylclabscm
RUN Rscript -e "install.packages('/build/ylclabscm', repos=NULL, type='source')"

# Wire exec scripts to PATH (R install.packages puts them in the library exec/ dir, not in PATH)
RUN exec_dir=$(Rscript -e "cat(system.file('exec', package='ylclabscm'))") && \
    ln -s "${exec_dir}/ylclabscm-batch-correct" /usr/local/bin/ylclabscm-batch-correct && \
    ln -s "${exec_dir}/ylclabscm-cluster" /usr/local/bin/ylclabscm-cluster && \
    ln -s "${exec_dir}/ylclabscm-diff-abundance" /usr/local/bin/ylclabscm-diff-abundance && \
    ln -s "${exec_dir}/ylclabscm-normalize" /usr/local/bin/ylclabscm-normalize && \
    ln -s "${exec_dir}/ylclabscm-reduce" /usr/local/bin/ylclabscm-reduce
