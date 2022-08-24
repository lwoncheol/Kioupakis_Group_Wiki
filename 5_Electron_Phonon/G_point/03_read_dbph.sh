#!/bin/bash

ftn -DMPI -o read_dbph read_dbph.f90

#for i in {1..106}
for i in {108..194}
do
    cd q$i/elph_dir
    cp ../../read_dbph .
    ./read_dbph
    cd ../../
done
