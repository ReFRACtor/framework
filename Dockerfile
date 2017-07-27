FROM centos:7

# Number of simultaneous build jobs for make
ENV BUILD_JOBS 10
ENV ABSCO_DIR /absco
ENV MERRA_DIR /merra
ENV BASE_DIR /refractor
ENV SRC_DIR $BASE_DIR/framework
ENV BUILD_DIR $BASE_DIR/build
ENV INSTALL_DIR $BASE_DIR/install

WORKDIR $BASE_DIR

# Get latest packages from Fedora Extra Packages for Enterprise Linux (EPEL) project
RUN yum install -y epel-release

# Install required packages for building
# Make a symbolic link for cmake3 since this is a non-standard name for the program
RUN yum install -y \
    file which \
    gcc gcc-gfortran gcc-c++ libtool rpm-build \
    cmake3 make patch zlib-devel bzip2-devel \
    hdf5 hdf5-devel readline-devel \
    python34 python34-devel python34-pip python34-numpy python34-nose \
    doxygen python-sphinx && \
    ln -s /usr/bin/cmake3 /usr/bin/cmake

# Install GSL 2+ from source since version in CentOS and EPEL are many years old in fact
ADD http://dl.fedoraproject.org/pub/fedora/linux/releases/26/Everything/source/tree/Packages/g/gsl-2.3-1.fc26.src.rpm /tmp/

RUN rpmbuild --rebuild /tmp/gsl-2.3-1.fc26.src.rpm && \
    yum install -y /root/rpmbuild/RPMS/x86_64/*.rpm && \
    rm /tmp/gsl-2.3-1.fc26.src.rpm

# Install Boost from source since we need >= 1.59
ADD https://dl.bintray.com/boostorg/release/1.64.0/source/boost_1_64_0.tar.gz /tmp/
RUN cd /tmp && tar zfvx boost_1_64_0.tar.gz && cd boost_1_64_0 && \
    ./bootstrap.sh --with-libraries=date_time,regex,iostreams,filesystem && \
    ./bjam install threading=single link=static,shared && \
    rm /tmp/boost_1_64_0.tar.gz

# Copy source code to container
ADD . $SRC_DIR

# Compile library and main executable
RUN mkdir -p $BUILD_DIR && \
    cd $BUILD_DIR && \
    cmake $SRC_DIR -DABSCO_DIR=$ABSCO_DIR -DMERRA_DIR=$MERRA_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR

RUN cd $BUILD_DIR && \
    make thirdparty 

RUN cd $BUILD_DIR && \
    make all -j $BUILD_JOBS && \
    make install
