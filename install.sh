#!/bin/bash

# Load environment variables from .env
source .env

mkdir -p build
pushd build
cmake -U LOG_LEVEL .. || exit 1
make || exit 1
popd