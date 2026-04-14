# yaffa
Yet Another Framework for Femtoscopy Analyses

## Installation:
Set the environment variables: `cp .env.example .env` and modify `.env` according to your machine.

1. Install the dependencies:
* `pip3 install [--user] pytest, dotenv, rich, numexpr, pylint`
* Install docker
* Install CECA
* Install yaml-cpp-devel
* Install python 3.10

1. OBSOLETE: add ` export YAFFA=/path/to/yaffa` to your .bashrc
1. run `bash install.sh`

1. OBSOLETE: If you use root inside the AliPhysics environment, add to the alienv script the line `export YAFFA=/path/to/yaffa`.
You can find the path to your alienv script with `which alienv`. If you use root outside of AliPhysics add the `YAFFA`
variable to your `.bashrc` file
