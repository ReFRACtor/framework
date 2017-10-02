FROM centos:7

# Number of simultaneous build jobs for make
ENV BUILD_JOBS 10
ENV ABSCO_DIR /data/absco
ENV MERRA_DIR /data/merra
ENV BASE_DIR /refractor
ENV SRC_DIR $BASE_DIR/framework
ENV BUILD_DIR $BASE_DIR/build
ENV INSTALL_DIR $BASE_DIR/install

WORKDIR $BASE_DIR

# Get latest packages from Fedora Extra Packages for Enterprise Linux (EPEL) project
RUN yum install -y epel-release && yum update -y

# Install required packages for building
# Make a symbolic link for cmake3 since this is a non-standard name for the program
RUN yum install -y \
    file which \
    gcc gcc-gfortran gcc-c++ libtool \
    cmake3 make patch zlib-devel bzip2-devel \
    hdf5 hdf5-devel readline-devel \
    python34 python34-devel python34-pip python34-numpy python34-nose \
    doxygen python-sphinx pcre-devel && \
    ln -s /usr/bin/cmake3 /usr/bin/cmake

# Install GSL 2+ from source since version in CentOS and EPEL are many years old in fact
ADD http://mirrors.ibiblio.org/gnu/ftp/gnu/gsl/gsl-2.4.tar.gz /tmp

RUN cd /tmp && tar zfvx gsl-2.4.tar.gz && cd gsl-2.4 && \
    ./configure && make && make install && \
    rm -rf /tmp/gsl-2.4*

# Install Boost from source since we need >= 1.59
ADD https://dl.bintray.com/boostorg/release/1.64.0/source/boost_1_64_0.tar.gz /tmp/
RUN cd /tmp && tar zfvx boost_1_64_0.tar.gz && cd boost_1_64_0 && \
    ./bootstrap.sh --with-libraries=date_time,regex,iostreams,filesystem && \
    ./bjam install threading=single link=static,shared && \
    rm -rf /tmp/boost_1_64_0*

# Install SWIG 3.0 from source, 2.0 is insufficient
ADD https://downloads.sourceforge.net/project/swig/swig/swig-3.0.12/swig-3.0.12.tar.gz /tmp
RUN cd /tmp && tar zfvx swig-3.0.12.tar.gz && cd swig-3.0.12 && \
    ./configure && make && make install && \
    rm -rf /tmp/swig-3.0.12*

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
