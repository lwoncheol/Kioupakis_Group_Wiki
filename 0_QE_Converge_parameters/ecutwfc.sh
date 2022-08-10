#!/bin/bash

for i in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100
do
    mkdir ecutwfc_$i
    cp batchscript ecutwfc_$i
    cp *UPF ecutwfc_$i
    cd ecutwfc_$i
    cat > scf.in <<EOF
&control
   prefix = 'BN'
   calculation = 'scf'
   restart_mode = 'from_scratch'
   wf_collect = .true.
   tstress = .true.
   tprnfor = .true.
   outdir = './'
   wfcdir = './'
   pseudo_dir = './'
   verbosity = 'high'
/
&system
   ibrav = 4
   celldm(1) = 4.731873875
   celldm(3) = 2.66022364217
   nat = 4
   ntyp = 2
   nbnd = 10
   ecutwfc = $i
/
&electrons
   electron_maxstep = 500
   conv_thr = 1.0d-10
   mixing_mode = 'plain'
   mixing_beta = 0.7
   mixing_ndim = 8
   diagonalization = 'david'
   diago_david_ndim = 4
   diago_full_acc = .true.
/
&ions
/
&cell
/
ATOMIC_SPECIES
   B  10.811   b.cpi.UPF
   N  14.007   n.cpi.UPF

ATOMIC_POSITIONS (crystal)
B       0.333333333   0.666666667   0.250000000
B       0.666666667   0.333333333   0.750000000
N       0.333333333   0.666666667   0.750000000
N       0.666666667   0.333333333   0.250000000

K_POINTS automatic
3 3 1 1 1 1
EOF
   sbatch batchscript
   cd ..
done
