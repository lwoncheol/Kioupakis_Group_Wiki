#!/bin/bash

module load berkeleygw/2.1.knl
module load python/2.7-anaconda-2019.07

echo -n "epsmat_hdf5_merge.py" > run_epsmat_merge.sh

for i in {1..76}
do
    echo -n " q$i/epsmat.h5" >> run_epsmat_merge.sh
done
    
    
