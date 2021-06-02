#!/bin/bash -l

# This serial script links the WFN files from the MF to the BGW folder.

# First, link all CD files
echo "Linking pseudopotential and charge densities"
for directory in {2,3,4,5,6}*-*/; do
  echo " Working on directory $directory"
  cd $directory
  mkdir GaN.save
  cd GaN.save
  cp ../../1-scf/GaN.save/charge-density.* .
  cp ../../1-scf/GaN.save/data-file-schema.xml .
  cd ..
  cp ../1-scf/*.UPF .
  cd ../
done
echo

echo "Linking WFNs to BGW directory"
cd ../2-bgw/

echo " Working on directory 1-epsilon"
cd 1-eps/
ln -sf ../../1-mf/2-wfn/wfn.cplx WFN
ln -sf ../../1-mf/3-wfnq/wfn.cplx WFNq
cd ../

echo " Working on directory 2-sigma"
cd 2-sig/
ln -sf ../../1-mf/4-wfn_co/rho.cplx ./RHO
ln -sf ../../1-mf/4-wfn_co/wfn.cplx WFN_inner
ln -sf ../../1-mf/4-wfn_co/vxc.dat .
ln -sf ../1-eps/eps0mat.h5 .
ln -sf ../1-eps/epsmat.h5 .
cd ..

echo " Working on directory 3-kernel"
cd 3-ker/
ln -sf ../../1-mf/4-wfn_co/wfn.cplx WFN_co
ln -sf ../1-eps/epsmat.h5 .
ln -sf ../1-eps/eps0mat.h5 .
cd ../

echo " Working on directory 4-absorption"
cd 4-abs/
ln -sf ../../1-mf/4-wfn_co/wfn.cplx WFN_co
ln -sf ../../1-mf/5-wfn_fi/wfn.cplx WFN_fi
ln -sf ../../1-mf/6-wfnq_fi/wfn.cplx WFNq_fi
ln -sf ../1-eps/epsmat.h5 .
ln -sf ../1-eps/eps0mat.h5 .
ln -sf ../2-sig/eqp1.dat eqp_co.dat
ln -sf ../3-ker/bsemat.h5 .
cd ../

echo
echo "Done!"
