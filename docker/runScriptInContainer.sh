#!/bin/bash

#   Attempts to execute python script specified by the argument within container.
#    The script path must be relative to current working directory, which will be
#    bind-mounted at /opt/xl/nadia-coecxs/NADIA/interfaces/ext within the container. 

echo target="$(pwd)"

docker run  -ti --rm  -v 'pwd':/opt/xl/nadia-coecxs/NADIA/interfaces/python/ext \
       -w /opt/xl/nadia-coecxs/NADIA/interfaces/python/ext \
       nadia_centos7:latest \
       "python $(1)" 


