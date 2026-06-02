#!/bin/bash
set -e -u

# Build the installer for MacOS.
# This script must be run from the root of the repository.
# Prerequisites: wheel has been build, e.g. using build_wheel.sh

WHL_NAME=$(cd dist && ls ./*.whl && cd ..)
pip install "dist/${WHL_NAME}[legacy]"

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/alphaviz.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
