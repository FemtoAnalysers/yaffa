#!/bin/bash

# Check if environment is properly set
if [[ ! -f .env ]]; then
    echo -e "\033[31mError: .env file not found! Exit!\033[0m"
    exit 1
fi
source .env

# Check that the path to YAFFA is defined
if [[ $YAFFA == '' ]]; then
    echo -e "\033[31mError: path to yaffa not found. Check your .env file. Exit!\033[0m"
    exit 1
fi
echo $YAFFA

if [[ $ALIENVLVL == '' ]]; then
    echo -e "\033[31mError: this script must be run inside the ALICE environment. Exit!\033[0m"
    exit 1
fi
