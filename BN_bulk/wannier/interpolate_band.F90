PROGRAM main
  IMPLICIT NONE
#ifdef MPI
  INCLUDE 'mpif.h'
#endif
  include 'fftw3.f'
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,256)

  !-- Edit before using:
  !-- Number of wannier, Number of bands, number of k-points (from wannier.win file)
  INTEGER, PARAMETER :: NW = 16, NB=40, NK1=12, NK2=12, NK3=6, NKTOT=NK1*NK2*NK3
  
  !-- Number of k-points ot interpolate
  INTEGER, PARAMETER :: nk = 1603!-- number of k-points to interpolate

  REAL(r8) ::  Elda(NB, NK1,NK2,NK3)
  COMPLEX(r8) :: Umat_local(Nb,NW,NK1,NK2,NK3), Umat(Nb,NW,NK1,NK2,NK3)
  COMPLEX(r8) :: u_matrix(NW,NW,NK1,NK2,NK3), u_matrix_opt(NB,NW,NK1,NK2,NK3)
  INTEGER :: ndimwin(NK1,NK2,NK3)
  COMPLEX(r8) :: H0(NB,NB,2,2,NK1,NK2,NK3), Hrot_local(NW,NW,2,2,NK1,NK2,NK3), Hrot(NW,NW,2,2,NK1,NK2,NK3),HR(NW,NW,2,2,NK1,NK2,NK3),Hint(NW,NW,2,2),Hint2(2*NW,2*NW)
	

  REAL(r8) :: Elda_out(2*NW,nk),Elda_out_local(2*NW,nk)
  REAL(r8) :: Kout(3),kx,ky,kz, kvbm(3),kcbm(3),evbm,ecbm,temp_x, temp_y, Rtemp(3),xco(nk),yco(nk),zco(nk)
  REAL(r8) :: E_VBM, E_CBM, band_x
  REAL(r8), PARAMETER :: ryd = 13.6056980659_r8, twopi =6.2831853071795864769_r8
  INTEGER :: map_local(3,125,NK1,NK2,NK3), map(3,125,NK1,NK2,NK3), ndegen_local(NK1,NK2,NK3), ndegen(NK1,NK2,NK3),nint
  REAL(r8) :: dmin_local(NK1,NK2,NK3), dmin(NK1,NK2,NK3), dtemp
  
  !-- FFTw:
  COMPLEX(r8) :: fftw_in( NK1,NK2,NK3),fftw_out(NK1,NK2,NK3)
  INTEGER(8) :: plan

  REAL(r8) :: A(3,3), tempkx, tempky, tempkz, dotprod

  !-- LAPACK
  INTEGER, PARAMETER :: LWORK=4*NW-1 !What doest 4 here stand for?
  INTEGER :: INFO
  complex(r8) WORK(LWORK)
  real(r8) :: RWORK(14*NW-2), W(2*NW) !What does 14 here stand for?

  REAL(r8) :: vtemp(3)
  INTEGER i,j,k,k1,k2,k3,m,n,idummy,is,isp,in,inp,isn,isnp,ix,iy,iz,kt1,kt2,kt3,kt(3),i1,i2,i3,ik,i_para
  INTEGER :: ixp1, ixm1, iyp1, iym1, izp1, izm1

  INTEGER :: ident, nprocs, ierr

#ifdef MPI
  INTEGER status(MPI_STATUS_SIZE)

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ident, ierr)
#else
  nprocs = 1
  ident = 0
#endif
  
  A(:,1) = (/  4.731873875000000e0, 0.000000000000000e0, 0.000000000000000e0 /)
  A(:,2) = (/ -2.365936937500000e0, 4.097922983254201e0, 0.000000000000000e0 /)
  A(:,3) = (/  0.000000000000000e0, 0.000000000000000e0, 1.258784275404157e1 /)

  !-- First pass: Find equivalent R-point in 1st BZ:
DO i_para = ident + 1, NK1 * NK2 * NK3, nprocs
  k1 = (i_para - 1) / (NK2 * NK3) + 1
  k2 = MOD((i_para - 1) / NK3, NK2) + 1
  k3 = MOD(i_para - 1, NK3) + 1
     vtemp = MATMUL( A, (/k1-1,k2-1,k3-1/))
     !dmin(k1,k2,k3) = SQRT(((k1-1)*A(1,1)+(k2-1)*A(1,2)+(k3-1)*A(1,3))**2 &
     !                   & +((k1-1)*A(2,1)+(k2-1)*A(2,2)+(k3-1)*A(2,3))**2 & 
     !                   & +((k1-1)*A(3,1)+(k2-1)*A(3,2)+(k3-1)*A(3,3))**2)
     dmin_local(k1,k2,k3) = SQRT(DOT_PRODUCT(vtemp,vtemp))
     map_local(:,1,k1,k2,k3) = (/ k1-1, k2-1, k3-1 /)
     DO ix = -2, 2; DO iy = -2, 2; DO iz = -2, 2
        kt = (/ k1+ix*NK1-1, k2+iy*NK2-1, k3+iz*NK3-1 /)
        !kt1 = k1 + ix * NK1 -1 
        !kt2 = k2 + iy * NK2 -1 
        !kt3 = k3 + iz * NK3 -1
            vtemp = MATMUL( A, kt)  
        dtemp = SQRT(DOT_PRODUCT(vtemp,vtemp))
        !dtemp = SQRT( (kt1*A(1,1)+kt2*A(1,2)+kt3*A(1,3))**2 &
        !     & +(kt1*A(2,1)+kt2*A(2,2)+kt3*A(2,3))**2 &
        !     & +(kt1*A(3,1)+kt2*A(3,2)+kt3*A(3,3))**2 )
        IF ( dtemp .LT. dmin_local(k1,k2,k3) ) THEN
           dmin_local(k1,k2,k3) = dtemp
           map_local(:,1,k1,k2,k3) = kt
        END IF
     END DO;     END DO; END DO
  END DO

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE( dmin_local, dmin, NK1 * NK2 * NK3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  dmin = dmin_local
#endif

  !-- Second pass: Find degeneracies

DO i_para = ident + 1, NK1 * NK2 * NK3, nprocs
  k1 = (i_para - 1) / (NK2 * NK3) + 1
  k2 = MOD((i_para - 1) / NK3, NK2) + 1
  k3 = MOD(i_para - 1, NK3) + 1
     ndegen_local(k1,k2,k3) = 0
     DO ix = -2, 2;     DO iy = -2, 2;     DO iz = -2, 2
        kt = (/ k1 + ix * NK1 -1 ,k2 + iy * NK2 -1 ,k3 + iz * NK3 -1 /)
            vtemp = MATMUL( A, kt)  
        dtemp = SQRT(DOT_PRODUCT(vtemp,vtemp))
        IF ( ABS(dtemp-dmin(k1,k2,k3)) .LT. 1.0E-6 ) THEN
           ndegen_local(k1,k2,k3) = ndegen_local(k1,k2,k3)+1
           map_local(:,ndegen_local(k1,k2,k3),k1,k2,k3) = kt
        END IF
     END DO;     END DO;     END DO
  END DO

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE( map_local, map, 3 * 125 * NK1 * NK2 * NK3, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE( ndegen_local, ndegen, NK1 * NK2 * NK3, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  map = map_local
  ndegen = ndegen_local
#endif

  OPEN(UNIT=1,FILE='u_matrix_opt.dat',ACTION='READ')
  READ(1,*)
  DO k1 = 1, NK1;  DO k2 = 1, NK2;  DO k3 = 1, NK3
     READ(1,*) ndimwin(k1,k2,k3)
     DO m = 1, NW
     DO n = 1, ndimwin(k1,k2,k3)
        READ(1,*) idummy,idummy,idummy,temp_x,temp_y
        u_matrix_opt(n,m,k1,k2,k3) = CMPLX(temp_x,temp_y)
     END DO
     END DO
  END DO;  END DO;  END DO
  CLOSE(1)

  OPEN(UNIT=1,FILE='u_matrix.dat',ACTION='READ')
  READ(1,*)
  DO k1 = 1, NK1;  DO k2 = 1, NK2;  DO k3 = 1, NK3
     DO m = 1, NW
     DO n = 1, NW
        READ(1,*) idummy,idummy,idummy,temp_x,temp_y
        u_matrix(n,m,k1,k2,k3) = CMPLX(temp_x,temp_y)
     END DO
     END DO
  END DO;  END DO;  END DO
  CLOSE(1)

  !-- Multiply matrices
DO i_para = ident + 1, NK1 * NK2 * NK3, nprocs
  k1 = (i_para - 1) / (NK2 * NK3) + 1
  k2 = MOD((i_para - 1) / NK3, NK2) + 1
  k3 = MOD(i_para - 1, NK3) + 1
     Umat_local(:,:,k1,k2,k3) = MATMUL(u_matrix_opt(:,:,k1,k2,k3),u_matrix(:,:,k1,k2,k3))
  END DO

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(Umat_local, Umat, Nb*NW*NK1*NK2*NK3, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  Umat = Umat_local
#endif

  !-- fiel with eigenvalues: GW or LDA
  open(unit=1,file='BN_GW.eig')

  !open(unit=1,file='sig.log.fullBZ')
  !open(unit=1,file='band.dat.out')
  DO k1 = 1, NK1;  DO k2 = 1, NK2;  DO k3 = 1, NK3
     DO in = 1, NB
        READ(1,*) idummy, idummy, elda(in,k1,k2,k3)
     END DO
  END DO;  END DO; END DO
  CLOSE(1)

!  OPEN(unit=1,file='somatrix.dat',ACTION='READ') 
  DO k1 = 1, NK1;  DO k2 = 1, NK2;  DO k3 = 1, NK3
!     READ(1,*)
!     READ(1,*)
     DO is  = 1, 2;     DO isp = 1, 2
        DO in  = 1, NB;        DO inp = 1, NB
!           read(1,*) idummy,idummy,idummy,idummy,temp_x,temp_y
!           H0(in,inp,is,isp,k1,k2,k3) = CMPLX(temp_x,temp_y)*ryd
           H0(in,inp,is,isp,k1,k2,k3) = CMPLX(0.0D0,0.0D0)
           IF ( is .EQ. isp .AND. in .EQ. inp ) THEN
              H0(in,inp,is,isp,k1,k2,k3) = H0(in,inp,is,isp,k1,k2,k3) + Elda(in,k1,k2,k3)
           END IF
        END DO;        END DO
     END DO;     END DO
  END DO;  END DO;  END DO
!  CLOSE(1)

  !-- rotate:
  PRINT*, 'Rotating...'
DO i_para = ident + 1, NK1 * NK2 * NK3, nprocs
  k1 = (i_para - 1) / (NK2 * NK3) + 1
  k2 = MOD((i_para - 1) / NK3, NK2) + 1
  k3 = MOD(i_para - 1, NK3) + 1
     DO is  = 1, 2;    DO isp = 1, 2
        DO in  = 1, NW; DO inp = 1, NW
           Hrot_local(in,inp,is,isp,k1,k2,k3) = CMPLX(0.0_r8,0.0_r8)
           DO n = 1, ndimwin(k1,k2,k3)
           DO m = 1, ndimwin(k1,k2,k3)
              Hrot_local(in,inp,is,isp,k1,k2,k3) = Hrot_local(in,inp,is,isp,k1,k2,k3) + &
                   & CONJG(Umat(n,in,k1,k2,k3)) * H0(n,m,is,isp,k1,k2,k3) * Umat(m,inp,k1,k2,k3)
           END DO
        END DO
     END DO; END DO
  END DO; END DO
END DO

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(Hrot_local, Hrot, NW*NW*2*2*NK1*NK2*NK3, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
 Hrot = Hrot_local
#endif
    
  PRINT*, 'Done rotating!'

  !-- FFT
  PRINT*, 'FFTing...'
  CALL dfftw_plan_dft_3d(plan,NK1,NK2,NK3,fftw_in,fftw_out, FFTW_FORWARD,FFTW_MEASURE)

  DO is  = 1, 2;  DO isp = 1, 2
     DO in  = 1, NW;     DO inp = 1, NW
        fftw_in(:,:,:) = Hrot(in,inp,is,isp,:,:,:)
        CALL dfftw_execute(plan)
        HR(in,inp,is,isp,:,:,:) = fftw_out(:,:,:)
     END DO;     END DO
  END DO;  END DO

  PRINT*, 'Done FFTing!'

  !-- interpolate in yz plane

  OPEN(UNIT=16,FILE='BN_band.kpt',ACTION='READ')
  READ(16,*)
  DO k = 1, nk
     READ(16,*) xco(k), yco(k), zco(k)
  END DO
  CLOSE(16)
  
  DO ik = ident+1, nk, nprocs

        DO is  = 1, 2
	DO isp = 1, 2
           DO in  = 1, NW
	   DO inp = 1, NW
              Hint(in,inp,is,isp) = CMPLX(0.0D0,0.0D0)
              DO k1 = 1, NK1
              DO k2 = 1, NK2
              DO k3 = 1, NK3
                 DO i = 1, ndegen(k1,k2,k3)
                    !-- Use R-vector in cartesian
                    !Rtemp(:) = map(1,i,k1,k2,k3)*A(:,1) + map(2,i,k1,k2,k3)*A(:,2) + map(3,i,k1,k2,k3)*A(:,3)
                    
                    dotprod = twopi*(xco(ik)*map(1,i,k1,k2,k3)+yco(ik)*map(2,i,k1,k2,k3)+zco(ik)*map(3,i,k1,k2,k3) )
                    !dotprod = xco(ik) * Rtemp(1) + yco(ik) * Rtemp(2) + zco(ik) * Rtemp(3)
                    temp_x = COS(dotprod)
                    temp_y = SIN(dotprod) !sign????
                    Hint(in,inp,is,isp) = Hint(in,inp,is,isp) + HR(in,inp,is,isp,k1,k2,k3) * &
                         & CMPLX(temp_x,temp_y)/(1.0_r8*NKTOT*ndegen(k1,k2,k3))
                 END DO
              END DO
              END DO
           END DO

        END DO
     END DO
  END DO
     END DO
	   
        !--diagonalize:

        !-- rearrange:
        DO is  = 1, 2
	DO isp = 1, 2
           DO in  = 1, NW
           DO inp = 1, NW
              isn  = in  + (is -1)*NW
              isnp = inp + (isp-1)*NW
              Hint2(isn,isnp) = Hint(in,inp,is,isp)
           END DO
           END DO
        END DO
        END DO

        CALL zheev('N','U',2*NW,Hint2,2*NW,W,WORK,LWORK,RWORK,INFO)

        DO n = 1, 2*NW
           Elda_out_local(n,ik) = W(n)
!           PRINT*, n, xco(ik), yco(ik), zco(ik), Elda_out_local(n,ik)
        END DO
     END DO

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE( Elda_out_local, Elda_out, 2*NW*nk, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  Elda_out = Elda_out_local
#endif

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_FINALIZE(ierr)
#endif

  OPEN(UNIT=17,FILE='BN_band.dat',ACTION='READ')
  OPEN(UNIT=18,FILE='out_band.dat',ACTION='WRITE')
  DO i = 1, NW
  DO k = 1, nk
     READ(17,*) band_x
     WRITE(18,*) band_x, (Elda_out(i * 2, k) + Elda_out(i * 2 - 1, k)) / 2
  END DO
!  READ(17,*) band_x
  WRITE(18,*) ""
  END DO
  CLOSE(17)
  CLOSE(18)

  STOP
END PROGRAM main
