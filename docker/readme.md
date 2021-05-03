### Instructions for using Nadia via Docker

To build the docker image:

1. from the directory containing the Dockerfile, run the command:

docker build -t nadia_centos7:latest .

2. Start an interactive bash shell in the container, mounting data and source code folders on the host:

docker run --name nadia -it \
       --mount type=bind,source="$(pwd)"/src,target=/opt/nadia-coecxs/NADIA/host/src \
       --mount type=bind,source="$(pwd)"/data,target=/opt/nadia-coecxs/NADIA/host/data \
       nadia_centos7:latest /bin/bash


docker run --name nadia -it --mount type=bind,source="$(pwd)"/src,target=/opt/nadia-coecxs/NADIA/host/src   --mount type=bind,source="$(pwd)"/data,target=/opt/nadia-coecxs/NADIA/host/data nadia_centos7:latest /bin/bash

3. Execute a command using exec, examples:

Run the planar CDI example:

docker exec -e LD_LIBRARY_PATH='/usr/local/hdf4/lib/:/usr/local/nadia/lib'  -it nadia bash -c 'python /opt/nadia-coecxs/NADIA/interfaces/python/planar_example.py'

Convert a file on the host using process in the container:

docker exec -e LD_LIBRARY_PATH='/usr/local/hdf4/lib/:/usr/local/nadia/lib'  -e PATH='/usr/local/nadia/bin' -it nadia '/usr/local/nadia/bin/tiff2ppm.exe host/data/object.tiff host/data/object.ppm'

### Compile and run code in host volume using container

Edit source files with .c extension in ./src on the host, which is bind mounted to /opt/nadia-coecxs/NADIA/host/src on the host.

Compile with 

   sh make.sh

 Note that it will generally be easier to use the interactive bash shell to compile and debug code.


### Utility Scripts

The following scripts can be used to simplify the above steps:

1. build.sh:   builds the docker image
2. run.sh:     runs the docker container
3. exec.sh <command>:  executes command within the container, e.g.:

    sh exec.sh '/usr/local/nadia/bin/tiff2ppm.exe host/data/object.tiff host/data/object.ppm'

