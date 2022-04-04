#!/bin/bash

echo "Copying and unzipping datasets to working directory"
cp /data/inputs/PreProcessedData.zip /code/PreProcessedData.zip
unzip /code/PreProcessedData.zip -d /code/PreProcessedData/
cp /data/inputs/UKCP18_subset.zip /code/UKCP18_subset.zip
unzip /code/UKCP18_subset.zip -d /code/UKCP18dir/

# echo "Running Hello World Model"
/usr/bin/mlrtapp/HEAT

# echo "Copying outputs to DAFNI outputs directory"
# cp /code/output.txt /data/outputs/output.txt

# echo "Printing contents of DAFNI outputs directory"
# ls -la /data/outputsd
