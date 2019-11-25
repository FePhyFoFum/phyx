#!/bin/bash

INPUTDIR=/home/josephwb/Work/Phylogenetics/phyx/src
TOOL=$INPUTDIR/pxs2fa


for f in *.phy; do
  echo -e "\nTesting file: " $f
  #more $f
  $TOOL -s $f
done
