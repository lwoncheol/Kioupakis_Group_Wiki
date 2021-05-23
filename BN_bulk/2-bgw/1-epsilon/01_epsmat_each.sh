#!/bin/bash

for i in {1..76}
do
    mkdir q$i
    cd q$i
    ln -sf ../../../1-mf/2-wfn/wfn.cplx WFN
    ln -sf ../../../1-mf/3-wfnq/wfn.cplx WFNq
    cp ../epsilon.inp .
    cp ../batchscript .
    awk "NR == $i" ../qpoints >> epsilon.inp
    echo 'end' >> epsilon.inp
    sbatch batchscript
    cd ..
done
