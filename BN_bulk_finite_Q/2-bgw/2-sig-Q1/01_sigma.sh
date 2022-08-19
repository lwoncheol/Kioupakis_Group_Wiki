#!/bin/bash

for i in {11..108}
do
    mkdir k$i
    cd k$i
    ln -sf ../../1-eps/q0/eps0mat.h5 .
    ln -sf ../../1-eps/epsmat.h5 .
    ln -sf ../../../1-mf/4-wfnq_co_Q1/wfn.cplx WFN_inner
    ln -sf ../../../1-mf/4-wfnq_co_Q1/rho.cplx RHO
    ln -sf ../../../1-mf/4-wfnq_co_Q1/vxc.dat .
    cp ../sigma.inp .
    cp ../batchscript .
    awk "NR == $((8*i-7))" ../kpoints >> sigma.inp
    awk "NR == $((8*i-6))" ../kpoints >> sigma.inp
    awk "NR == $((8*i-5))" ../kpoints >> sigma.inp
    awk "NR == $((8*i-4))" ../kpoints >> sigma.inp
    awk "NR == $((8*i-3))" ../kpoints >> sigma.inp
    awk "NR == $((8*i-2))" ../kpoints >> sigma.inp
    awk "NR == $((8*i-1))" ../kpoints >> sigma.inp
    awk "NR == $((8*i))" ../kpoints >> sigma.inp
    echo 'end' >> sigma.inp
    sbatch batchscript
    cd ..
done
