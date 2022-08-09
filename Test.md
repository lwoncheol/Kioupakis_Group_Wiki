# Projected Density Of States (PDOS)

<aside>
💡 These instructions cover how to calculate projected density of state (PDOS) using **“projwfc.x”** executable of Quantum Espresso

</aside>

# 0. Overview

- For non-relativistic pseudopotentials (spin-orbit coupling not considered), wavefunctions are projected on to the atomic orbitals (s, p, d, and so on).
- For fully-relativistic pseudopotentials (spin-orbit coupling included), spinor wavefunctions are projected on to the eigenstates of total-angular momentum (j)
- (Special case) If you want to project spinors on to the atomic orbitals, check section 3.
- For the detailed description for each flag, check: [https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html)

# 1. N**on-relativistic pseudopotential (w/o spin-orbit)**

1. Run scf calculation
    
    ```fortran
    &control
       prefix = 'PbI2'
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
       celldm(1) = 8.613371056
       celldm(3) = 1.5326897762
       nat = 3
       ntyp = 2
       nbnd = 24
       ecutwfc = 90.0
    /
    &electrons
       electron_maxstep = 500
       conv_thr = 1.0d-10
       mixing_mode = 'plain'
       mixing_beta = 0.5
       mixing_ndim = 8
       diagonalization = 'david'
       diago_david_ndim = 2
       diago_full_acc = .true.
    /
    &ions
    /
    &cell
    /
    ATOMIC_SPECIES
       Pb  207.2      Pb.cpi.UPF
       I   126.90447   I.cpi.UPF
    
    ATOMIC_POSITIONS (crystal)
    Pb       0.000000000   0.000000000   0.000000000
    I        0.333333333   0.666666667   0.267500000
    I        0.666666667   0.333333333   0.732500000
    
    K_POINTS automatic
    12 12 8 1 1 1
    ```
    
2. Run nscf (or bands) calculations. Include bands and kpoints you want to analyze 
    
    ```fortran
    &control
       prefix = 'PbI2'
       calculation = 'bands'
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
       celldm(1) = 8.613371056
       celldm(3) = 1.5326897762
       nat = 3
       ntyp = 2
       nbnd = 24
       ecutwfc = 90.0
    /
    &electrons
       electron_maxstep = 500
       conv_thr = 1.0d-10
       mixing_mode = 'plain'
       mixing_beta = 0.5
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
       Pb  207.2      Pb.cpi.UPF
       I   126.90447   I.cpi.UPF
    
    ATOMIC_POSITIONS (crystal)
    Pb       0.000000000   0.000000000   0.000000000
    I        0.333333333   0.666666667   0.267500000
    I        0.666666667   0.333333333   0.732500000
    
    K_POINTS crystal_b
    8
     0.0000000000 0.0000000000 0.0000000000  200
     0.3333333333 0.3333333333 0.0000000000  98
     0.3333333333 0.3333333333 0.5000000000  200
     0.0000000000 0.0000000000 0.5000000000  98
     0.0000000000 0.0000000000 0.0000000000  173
     0.5000000000 0.0000000000 0.0000000000  98
     0.5000000000 0.0000000000 0.5000000000  173
     0.0000000000 0.0000000000 0.5000000000  1
    ```
    
3. (You can skip this step) Run `bands.x` executable (implemented in Quantum Espresso) to print out the band structure
    
    ```fortran
    &bands
       prefix  = 'PbI2'
       outdir = './'
       filband = 'bands.dat'
    /
    ```
    
4. Run `projwfc.x` executable (implemented in Quantum Espresso)
    
    ```fortran
    &projwfc
        outdir='./'
        prefix='PbI2'
        Emin = 0.0, Emax = 10.0, DeltaE = 0.01 % energy grid step in eV
        ngauss=0, degauss=0.001                % gaussian broadening, Ry (not eV)
        kresolveddos = .true.
    /
    ```
    
    - **DeltaE** and **degauss** are the most important parameters. They decide the precision and the resolution of PDOS calculation. The values specified above guarantee high-precision, but the size of the output files will be huge. Do not increase these values unless you are doing extremely heavy calculations.
5. Check the output files. For example, followings are output files for bulk PbI2. There are total 3 atoms (1 Pb and 2 I atoms) in a unit cell, and each atom has s,p,d orbitals.
    
    ```fortran
    % Output script file
    projwfc.out
    
    % PDOS for each orbital of each atom
    PbI2.pdos_atm#1(Pb)_wfc#1(s)       
    PbI2.pdos_atm#1(Pb)_wfc#2(p)
    PbI2.pdos_atm#1(Pb)_wfc#3(d)
    PbI2.pdos_atm#2(I)_wfc#1(s)
    PbI2.pdos_atm#2(I)_wfc#2(p)
    PbI2.pdos_atm#2(I)_wfc#3(d)
    PbI2.pdos_atm#3(I)_wfc#1(s)
    PbI2.pdos_atm#3(I)_wfc#2(p)
    PbI2.pdos_atm#3(I)_wfc#3(d)
    
    % Total PDOS
    PbI2.pdos_tot
    ```
    
6. The format of output files is:
    
    
    | #ik | E (eV) | LDOS (E) | PDOS (E) | PDOS(E) | … |
    | --- | --- | --- | --- | --- | --- |
    | 1 | 0.000 | 0.000E+00 | 0.000E+00 | 0.000E+00 | … |
    | 1 | 0.010 | 0.000E+00 | 0.000E+00 | 0.000E+00 | … |
7. Orbital orders in each output file follows:
    
    ```fortran
    for l=1:
      1 pz     (m=0)
      2 px     (real combination of m=+/-1 with cosine)
      3 py     (real combination of m=+/-1 with sine)
    
    for l=2:
      1 dz2    (m=0)
      2 dzx    (real combination of m=+/-1 with cosine)
      3 dzy    (real combination of m=+/-1 with sine)
      4 dx2-y2 (real combination of m=+/-2 with cosine)
      5 dxy    (real combination of m=+/-2 with sine)
    
    ...
    ```
    

# 2. Fully**-relativistic pseudopotential (w/ spin-orbit)**

1. Run scf calculation. Check the use of `noncolin = .true.` ,`constrained_magnetization = ‘none’`, `lsipnorb = .true.`
    
    ```fortran
    &control
       prefix = 'PbI2'
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
       celldm(1) = 8.613371056
       celldm(3) = 1.5326897762
       nat = 3
       ntyp = 2
       ecutwfc = 90.0
       input_dft = 'pz'
       **noncolin = .true.
       constrained_magnetization='none'
       lspinorb = .true.**
    /
    &electrons
       electron_maxstep = 500
       conv_thr = 1.0d-10
       mixing_mode = 'plain'
       mixing_beta = 0.5
       mixing_ndim = 8
       diagonalization = 'david'
       diago_david_ndim = 2
       diago_full_acc = .true.
    /
    &ions
    /
    &cell
    /
    ATOMIC_SPECIES
       Pb  207.2      Pb.UPF
       I   126.90447   I.UPF
    
    ATOMIC_POSITIONS (crystal)
    Pb       0.000000000   0.000000000   0.000000000
    I        0.333333333   0.666666667   0.267500000
    I        0.666666667   0.333333333   0.732500000
    
    K_POINTS automatic
    12 12 8 1 1 1
    ```
    
2. Run nscf (or bands) calculations. Include bands and kpoints you want to analyze 
    
    ```fortran
    &control
       prefix = 'PbI2'
       calculation = 'bands'
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
       celldm(1) = 8.613371056
       celldm(3) = 1.5326897762
       nat = 3
       ntyp = 2
       nbnd = 48
       ecutwfc = 90.0
       input_dft = 'pz'
       **noncolin = .true.
       constrained_magnetization='none'
       lspinorb = .true.**
    /
    &electrons
       electron_maxstep = 500
       conv_thr = 1.0d-10
       mixing_mode = 'plain'
       mixing_beta = 0.5
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
       Pb  207.2      Pb.UPF
       I   126.90447   I.UPF
    
    ATOMIC_POSITIONS (crystal)
    Pb       0.000000000   0.000000000   0.000000000
    I        0.333333333   0.666666667   0.267500000
    I        0.666666667   0.333333333   0.732500000
    
    K_POINTS crystal_b
    8
     0.0000000000 0.0000000000 0.0000000000  200
     0.3333333333 0.3333333333 0.0000000000  98
     0.3333333333 0.3333333333 0.5000000000  200
     0.0000000000 0.0000000000 0.5000000000  98
     0.0000000000 0.0000000000 0.0000000000  173
     0.5000000000 0.0000000000 0.0000000000  98
     0.5000000000 0.0000000000 0.5000000000  173
     0.0000000000 0.0000000000 0.5000000000  1
    ```
    
3. (You can skip this step) Run `bands.x` executable (implemented in Quantum Espresso) to print out the band structure
    
    ```fortran
    &bands
       prefix  = 'PbI2'
       outdir = './'
       filband = 'bands.dat'
    /
    ```
    
4. Run `projwfc.x` executable (implemented in Quantum Espresso)
    
    ```fortran
    &projwfc
        outdir='./'
        prefix='PbI2'
        Emin = 0.0, Emax = 10.0, DeltaE = 0.01 % energy grid step in eV
        ngauss=0, degauss=0.001                % gaussian broadening, Ry (not eV)
        kresolveddos = .true.
    /
    ```
    
    - **DeltaE** and **degauss** are the most important parameters. They decide the precision and the resolution of PDOS calculation. The values specified above guarantee high-precision, but the size of the output files will be huge. Do not increase these values unless you are doing extremely heavy calculations.
5. Check the output files. For example, followings are output files for bulk PbI2. There are total 3 atoms (1 Pb and 2 I atoms) in a unit cell, and each atom has s,p,d orbitals.
    
    ```fortran
    % Output script file
    projwfc.out
    
    % PDOS for each orbital of each atom
    PbI2.pdos_atm#1(Pb)_wfc#1(d_j2.5) 
    PbI2.pdos_atm#1(Pb)_wfc#2(d_j1.5)
    PbI2.pdos_atm#1(Pb)_wfc#3(s_j0.5)
    PbI2.pdos_atm#1(Pb)_wfc#4(p_j1.5)
    PbI2.pdos_atm#1(Pb)_wfc#5(p_j0.5)
    PbI2.pdos_atm#2(I)_wfc#1(s_j0.5)
    PbI2.pdos_atm#2(I)_wfc#2(p_j1.5)
    PbI2.pdos_atm#2(I)_wfc#3(p_j0.5)
    PbI2.pdos_atm#2(I)_wfc#4(d_j2.5)
    PbI2.pdos_atm#2(I)_wfc#5(d_j1.5)
    PbI2.pdos_atm#3(I)_wfc#1(s_j0.5)
    PbI2.pdos_atm#3(I)_wfc#2(p_j1.5)
    PbI2.pdos_atm#3(I)_wfc#3(p_j0.5)
    PbI2.pdos_atm#3(I)_wfc#4(d_j2.5)
    PbI2.pdos_atm#3(I)_wfc#5(d_j1.5)
    
    % Total PDOS
    PbI2.pdos_tot
    ```
    
6. The format of output files is:
    
    
    | #ik | E (eV) | LDOS (E) | PDOS_1 (E) | PDOS_2(E) | … |
    | --- | --- | --- | --- | --- | --- |
    | 1 | 0.000 | 0.000E+00 | 0.000E+00 | 0.000E+00 | … |
    | 1 | 0.010 | 0.000E+00 | 0.000E+00 | 0.000E+00 | … |
7. Orbital orders in each output file follows: **(Need to double check)**
    
    ```fortran
    for j=0.5:
      1 mj = -0.5
      2 mj = +0.5
    
    for j=1.5:
      1 mj = -1.5
      2 mj = -0.5
      3 mj = +0.5
      4 mj = +1.5
    
    ...
    ```
    

# 3. Special case: How to project spinors on the atomic orbitals

- Quantum Espresso doesn’t support this type of calculation, but still you can cheat in projecting your eigenstates into the non-relatistic atomic states ([https://lists.quantum-espresso.org/pipermail/users/2021-July/047654.html](https://lists.quantum-espresso.org/pipermail/users/2021-July/047654.html))

1. Preceding steps (scf, bands, bands.x) are the same as Section 2.
2. **(Very Important)** Open the `data-file-schema.xml` file stored in the [`prefix.save`](http://prefix.save) folder. You will find the **magnetization** section. Change the **spinorbit** element from true to false.
    
    ```fortran
    <magnetization>
          <lsda>false</lsda>
          <noncolin>true</noncolin>
          <spinorbit>**true**</spinorbit>       % <-- Change this to false !
          <total>0.000000000000000e0</total>
          <absolute>0.000000000000000e0</absolute>
          <do_magnetization>false</do_magnetization>
    </magnetization>
    ```
    
3. Run `projwfc.x` executable. Following steps are the same as Section 2.
