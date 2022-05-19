#!/bin/bash

echo "Copying and unzipping datasets to working directory"
cp /data/inputs/PreProcessedData/PreProcessedData.tar.gz /code/PreProcessedData.tar.gz
tar -xf /code/PreProcessedData.tar.gz -C /code/PreProcessedData/
cp /data/inputs/UKCP18dir/UKCP18_subset.tar.gz /code/UKCP18_subset.tar.gz
tar -xf /code/UKCP18_subset.tar.gz -C /code/UKCP18dir/
cp /data/inputs/input_file/* /code/


# echo "Running containerised model"
/usr/bin/mlrtapp/HEAT
