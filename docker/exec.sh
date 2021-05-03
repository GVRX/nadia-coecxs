#!/bin/bash

docker exec -e LD_LIBRARY_PATH='/usr/local/hdf4/lib/:/usr/local/nadia/lib'  -it nadia $1 
