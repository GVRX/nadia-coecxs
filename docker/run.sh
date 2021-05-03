#!/bin/bash

docker run --name nadia -it --rm --mount type=bind,source="$(pwd)"/src,target=/opt/nadia-coecxs/NADIA/host/src   --mount type=bind,source="$(pwd)"/data,target=/opt/nadia-coecxs/NADIA/host/data nadia_centos7:latest /bin/bash


