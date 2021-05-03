FROM centos:7

LABEL maintainer Grant van Riessen <g.vanriessen@latrobe.edu.au>

#RUN curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo \
RUN yum -y  update \
    && yum -y groupinstall  "development tools" \
    && yum -y install epel-release \
    && yum -y install  zlib-devel \
                      bzip2-devel \
                      openssl-devel \
                      ncurses-devel \
                      sqlite-devel \
                      readline-devel \
                      tk-devel \
                      gdbm-devel \
                      db4-devel \
                      libpcap-devel \
                      xz-devel \
                      expat-devel \
                      wget \
                      openmpi-devel \
                      fftw3-devel \
                      libtiff-devel \
                      lapack-devel \
                      blas-devel \
                      mlocate \
                      libjpeg-devel \
                      libffi-devel \
                      mpich-devel \
                      libzip \
                      gperftools-devel \
                      gperftools \
                      Cython \
                      python-devel

ENV LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
#ENV CC=/usr/lib64/openmpi/bin/mpicc
ENV /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /tmp
COPY requirements.txt requirements.txt
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py \    
    && python get-pip.py \ 
    && rm get-pip.py \
    && pip install -r requirements.txt \
    && rm requirements.txt

WORKDIR /opt
RUN wget https://support.hdfgroup.org/ftp/HDF/releases/HDF4.2.14/src/hdf-4.2.14.tar.gz \
    && tar -zxvf hdf-4.2.14.tar.gz \
    && rm hdf-4.2.14.tar.gz \
    && cd hdf-4.2.14 \
    && ./configure --prefix=/usr/local/hdf4 --with-zlib --with-jpeg --enable-shared --disable-netcdf  --disable-fortran \
    && make -j$threads \
    && make install

WORKDIR /opt
ARG HDF=/usr/local/hdf4
ENV LD_LIBRARY_PATH=/opt/nadia-coecxs/NADIA/interfaces/python/pyNADIA:$HDF:$LD_LIBRARY_PATH

RUN git clone https://github.com/GVRX/nadia-coecxs.git\
    && cd nadia-coecxs/NADIA \
    && mkdir lib \
    && mkdir bin \
    && mkdir interfaces/python/ext \
    && export DO_IDL=FALSE \
    && ./configure  --prefix=/usr/local/nadia --with-mfhdf_lib=$HDF/lib/ --with-mfhdf_inc=$HDF/include \ 
    && make
