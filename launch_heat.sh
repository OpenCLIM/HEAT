#!/bin/bash

echo "Copying and unzipping datasets to working directory"
cp /data/inputs/PreProcessedData.zip /code/PreProcessedData.zip
cp /data/inputs/UKCP18_subset.zip /code/UKCP18_subset.zip

echo "Running Hello World Model"
/usr/bin/mlrtapp/HEAT

echo "Copying outputs to DAFNI outputs directory"
# cp /code/output.txt /data/outputs/output.txt

echo "Printing contents of DAFNI outputs directory"
# ls -la /data/outputsd
