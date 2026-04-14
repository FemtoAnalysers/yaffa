#!/bin/bash

docker run -e RUN_LOCAL=true -e LOG_LEVEL=WARN -e VALIDATE_PYTHON_PYLINT=true -v $PWD:/tmp/lint github/super-linter \
    &> >(grep -v "fatal: not a git repository (or any of the parent directories): .git")