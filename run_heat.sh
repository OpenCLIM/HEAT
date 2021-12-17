#!/bin/bash



echo "Running HEAT"
/usr/bin/mlrtapp/heatdocker


echo "Printing contents of DAFNI outputs directory"
ls -la /data/outputs
