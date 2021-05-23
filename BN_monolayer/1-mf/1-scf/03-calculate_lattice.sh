#!/bin/bash -l

module load berkeleygw

surface.x surface.inp &> lattice.out

# get x,y,z lengths from lattice.out
  export dx2=`grep dx lattice.out | awk '{print $9;}'`
  export dy2=`grep dy lattice.out | awk '{print $9;}'`
  export dz2=`grep dz lattice.out | awk '{print $9*2;}'`

  echo "How large does the c-axis have to be so that c/2 contains 99% of the charge density?"
  echo $dz2 "a.u."
