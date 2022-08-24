#!/bin/bash

cd global_phonon
for i in {1..194}
do
    mkdir q$i
    cd q$i
    cp ../../global_scf/*.UPF .
    mkdir Si.save
    cp -R ../../global_scf/Si.save .
    cat > ph.in <<EOF
phonons for Si
&inputph
 outdir = '.'
 prefix = 'Si'
 trans = .true.
 epsil = .false.
 ldisp = .false.
 qplot = .true.
 electron_phonon = 'dvscf'
 fildyn = 'Si.dyn'
 fildvscf = 'dvscf'
 tr2_ph = 1.0d-14
/
1
EOF
    awk "NR == $i" ../../qpoints.dat >> ph.in
    cat > batchscript <<EOF    
#!/bin/bash -l
#SBATCH -J Si-$i
#SBATCH -p regular
##SBATCH -q flex
##SBATCH --time-min=00:50:00
##SBATCH --time=48:00:00
#SBATCH --error=err_%j
#SBATCH --output=out_%j
#SBATCH -N 1
#SBATCH -S 4
#SBATCH -C knl,quad,cache
#SBATCH -t 00:49:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woncheol@umich.edu

srun -n 64 bash -c '/global/cfs/cdirs/m1380/modules/cori/q-e-6.7-auger-hacks/bin/ph.x -npool 64 < ph.in > ph.out'
EOF

    sbatch batchscript
    cd ..
done
cd ..
