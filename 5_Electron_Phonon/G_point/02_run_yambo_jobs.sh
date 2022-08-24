#!/bin/bash

for i in {1..194}
do
    mkdir q$i
    cd q$i
    
    cp -R ../global_nscf/Si.save/ .
    cp ../../global_phonon/q$i/Si.dyn1 .
    cp ../../global_phonon/q$i/*.UPF .
    cp -R ../../global_phonon/q$i/_ph0 .
    
    cat > ph.in <<EOF
phonons for Si
&inputph
 outdir = '.'
 prefix = 'Si'
 trans = .false.
 epsil = .false.
 ldisp = .false.
 qplot = .true.
 electron_phonon = 'yambo'
 fildyn = 'Si.dyn'
 fildvscf = 'dvscf'
 tr2_ph = 1.0d-14
/
1
EOF
    awk "NR == $i" ../../qpoints.dat >> ph.in
    cat > batchscript <<EOF    
#!/bin/bash -l
#SBATCH -J Si-yambo-$i
#SBATCH -p regular
##SBATCH -q flex
##SBATCH --time-min=00:50:00
##SBATCH --time=48:00:00
#SBATCH --error=err_%j
#SBATCH --output=out_%j
#SBATCH -N 1
#SBATCH -S 67
#SBATCH -C knl,quad,cache
#SBATCH -t 00:09:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woncheol@umich.edu

srun -n 1 bash -c '/global/cfs/cdirs/m1380/modules/cori/q-e-6.7-auger-hacks/bin/ph.x -npool 1 < ph.in > ph.out'
EOF

    sbatch batchscript
    cd ..
done
