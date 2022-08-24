!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    Read electron-phonon matrix elements in Bloch representation (*.epb, computed from EPW)
!!
!!    Reading code is from elphon_shuffle_wrap.f90
!!
!!    Woncheol Lee, 04/24/2022
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program read_dbph

  !#include "f_defs.h"

  implicit none


  integer, parameter :: DP = kind(1.0d0)

!!!INTEGER, PARAMETER :: qe_dp = selected_real_kind(14,200)

  character(len=100) :: fname
  !! string of input file name
  
  integer :: nbnd
  !! number of bands considered in nscf calculation

  integer :: nmodes
  !! number of phonon branches (3 * number of atoms in a unit cell)

  integer :: nksq
  !! number of kpoints considered in nscf calculation
  
  complex(DP), allocatable :: g(:,:,:,:) 
  !! electron-phonon matrix (nkpt, ibnd, jbnd, nmode)

  real(DP), allocatable :: omega(:)
  !! phonon frequency in unit of Ry

  real(DP), allocatable :: omega2(:)
  !! phonon frequency stored in db files. Unit is Ry**2
  
  logical :: err
  
  integer :: i, j, k, m
  !! loop integers
  
  integer :: iug = 125
  !! ID for e-ph matrix output file

  integer :: iuw = 135
  !! ID for phonon frequency output file

  integer ::iub = 145
  !! ID for band energy output file

  real(DP) :: qpt(3), fdummy
  real(DP), allocatable :: kpts_cart(:,:)
  real(DP), allocatable :: dummy_vec(:,:,:)
  real(DP), allocatable :: dummy_nbnd(:)

  real(DP) :: ryd2ha = 0.5
  real(DP) :: ryd2eV = 13.6056980659

!!!!!!!!!!!!!!!

  fname = 's.dbph_000001'
  nksq = 1
  nbnd = 8
  nmodes = 6
  
  allocate(g(nksq,nbnd,nbnd,nmodes))
  allocate(omega(nmodes))
  allocate(omega2(nmodes))
  allocate(kpts_cart(3,nksq))
  allocate(dummy_vec(nmodes,nmodes/3,3))
  allocate(dummy_nbnd(nbnd))

!!!!!!!!!!!!!!!
  
  open(unit=1, file=fname, action='READ',form="unformatted")
  !read header
  read(1) nmodes, nksq, nbnd
  
  !continue reading header
  read(1) fdummy,qpt,kpts_cart !fdummy is alat in Bohr
  !may want to check kpts_cart vs expected coordinate for sanity
  read(1) omega2(:) !Ry**2
  omega = sqrt(omega2) !Ry
  
  !done reading header
  !now read g (and other data as dummy) for each kpt
  do i = 1, nksq
     read(1) g(i,:,:,:) !Ry^(3/2) (still need to divide by sqrt(2*omega)
     do j = 1, nbnd
        do k = 1, nbnd
           g(i,j,k,:) = g(i,j,k,:)/sqrt(2*omega(:)) !convert to Ry, hbar is 1 in Ry too, not included here
        end do
     end do
     read(1) dummy_vec(:,:,:)
     read(1) dummy_nbnd(:)
     read(1) dummy_nbnd(:)
  end do
  !convert to Ha to match internal units of auger code
  !g = g*ryd2ha
  !omega = omega*ryd2ha

  g = g*ryd2eV
  omega = omega*ryd2eV
  err = .false.
  close(1)

  dummy_nbnd = dummy_nbnd*ryd2eV
!  PRINT*, dummy_nbnd

  OPEN(iug, FILE = 'g_matrix_G_VB_11.dat', ACTION = 'WRITE')
!  WRITE(iug, '(4a10, 3a15)') 'ik', 'ibnd', 'jbnd', 'imodes', 'real(g)', 'imag(g)', 'omega'

  DO i = 1, nksq
     DO j = 1, 1
     !DO j = 1, nbnd
        DO k = 1, 1
       !DO k = 1, nbnd    
           DO m = 1, nmodes
              WRITE(iug, '(1e15.6)') abs(g(i,j,k,m))*1000
              !WRITE(iug, '(4i10, 3e15.6)') i, j, k, m, real(g(i,j,k,m)), aimag(g(i,j,k,m)), omega(m)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  CLOSE(iug)

  
  OPEN(iuw, FILE = 'phfreq.dat', ACTION = 'WRITE')

  DO i = 1, nksq
     DO j = 1, 1
     !DO j = 1, nbnd
        DO k = 1, 1
       !DO k = 1, nbnd
           DO m = 1, nmodes
              WRITE(iuw, '(1e15.6)') omega(m)*1000
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  CLOSE(iuw)

  OPEN(iub, FILE = "bands.dat", ACTION = 'WRITE')

  DO i = 1, nbnd
     WRITE(iub, '(1e15.6)') dummy_nbnd(i)
  ENDDO

  CLOSE(iub)
  
end program read_dbph
