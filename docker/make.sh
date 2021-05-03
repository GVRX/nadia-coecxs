#!/bin/bash

docker exec -e LD_LIBRARY_PATH='/usr/local/hdf4/lib/:/usr/local/nadia/lib' -w /opt/nadia-coecxs/NADIA/host/src -it nadia make
