FROM continuumio/miniconda3

# Number of simultaneous build jobs for make
ENV BUILD_JOBS 5
ENV ABSCO_DIR /data/absco
ENV MERRA_DIR /data/merra
ENV BASE_DIR /refractor
ENV SRC_DIR $BASE_DIR/framework
ENV BUILD_DIR $BASE_DIR/build
ENV INSTALL_DIR $BASE_DIR/install

WORKDIR $BASE_DIR

# Make RUN commands use `bash --login` so conda works correctly
# Meaning that /root/.bashrc gets evaluated to initalize conda
SHELL ["/bin/bash", "--login", "-c"]

# Copy source code to container
ADD . $SRC_DIR

# Install mamba to avoid long environment solving times from conda
RUN conda install -y -c conda-forge mamba 

# Set up compilation tools and packages
RUN mamba env update -n base --file $SRC_DIR/environment.yml

# Compile library and main executable
# Need to refresh conda environment variables from install of environment
RUN cmake -S$SRC_DIR -B$BUILD_DIR -DABSCO_DIR=$ABSCO_DIR -DMERRA_DIR=$MERRA_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR

# Build ReFRACtor Framework
RUN cd $BUILD_DIR && \
    make -j $BUILD_JOBS install
