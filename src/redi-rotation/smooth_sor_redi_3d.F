#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      SUBROUTINE SMOOTH_SOR_REDI_3D (
     U     fld_in,
     I     smoothOpNb,
     I     center, east, west, north, south, up, down,
     I     east_up, east_down, west_up, west_down,
     I     north_up, north_down, south_up, south_down,
     I     smooth3Dmask,
     I     mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_jacobi_3d
C     | o Invert matrix A via jacobi iterations to solve for u,
C     |   given right hand side b:
C     |          Au = b
C     |
C     | o matrix A is "Laplacian like" and defined by
C     |   center, north, south, east, west fields
C     |
C     | o b := fld_in -> fld_rhs
C     |   u:= fld_out (fld_in is overwritten with this result at end)
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
C   center, north, south, east, west, up, down :
C       fields that define the laplacian stencil at each point
C       for a point (i,j,k)
C           center = (i,j,k)
C           east  = (i+1, j  ,  k  ); west  = (i-1, j  , k  )
C           north = (i  , j+1,  k  ); south = (i  , j-1, k  )
C           up    = (i  , j  ,  k+1); down  = (i  , j  , k-1)
C
      _RL fld_in      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER smoothOpNb
      _RL center      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL east_up   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west_up   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north_up  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south_up  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east_down (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west_down (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north_down(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south_down(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS smooth3Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER myThid

C --- Local parameters
      CHARACTER*(max_len_mbuf) msgBuf
      INTEGER i,j,k, bi, bj
      INTEGER jacobi_iter
      _RL TOL, TOL_SQ, fld_in_norm, relErr, absErr
      _RL sor

      _RL tile_in_norm(nSx,nSy)
      _RL tile_absErr(nSx,nSy)

      _RL fld_rhs     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL fld_out     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL recip_center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL Ax          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

C --- Initialize solver stuff
      jacobi_iter = 0
      TOL = smooth3DEllipticTol(smoothOpNb)
      TOL_SQ = TOL*TOL
      relErr = 100. _d 0
      fld_in_norm = zeroRL
      absErr = zeroRL
      sor = smooth3DSOROmega(smoothOpNb)

      do bj = myByLo(myThid),myByHi(myThid)
       do bi = myBxLo(myThid),myBxHi(myThid)
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           fld_rhs(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)
           fld_in (i,j,k,bi,bj) = oneRL
           fld_out(i,j,k,bi,bj) = oneRL
           if (ABS(center(i,j,k,bi,bj)).gt.zeroRL) then
            recip_center(i,j,k,bi,bj) = oneRL / center(i,j,k,bi,bj)
           else
            recip_center(i,j,k,bi,bj) = zeroRL
           endif
          enddo
         enddo
        enddo
       enddo
      enddo

C --- Precompute the error denominator once
      do bj = myByLo(myThid),myByHi(myThid)
       do bi = myBxLo(myThid),myBxHi(myThid)
        tile_in_norm(bi,bj) = zeroRL
        do k = 1,Nr
         do j = 1,sNy
          do i = 1,sNx
           tile_in_norm(bi,bj)  = tile_in_norm(bi,bj) +
     &          fld_rhs(i,j,k,bi,bj)**2
          enddo
         enddo
        enddo
       enddo
      enddo

      CALL GLOBAL_SUM_TILE_RL(tile_in_norm, fld_in_norm, myThid)

C --- The jacobi iterations
      do while ((jacobi_iter.lt.smooth3DJacobiMaxIters(smoothOpNb))
     &     .and.(relErr.ge.TOL_SQ))

       ! reset error and update counter
       relErr      = zeroRL
       jacobi_iter = jacobi_iter + 1

       do bj = myByLo(myThid),myByHi(myThid)
        do bi = myBxLo(myThid),myBxHi(myThid)
         tile_absErr (bi,bj) = zeroRL
         do k = 1,Nr
          do j = 1,sNy
           do i = 1,sNx

            ! Update "final" estimate
            fld_out(i,j,k,bi,bj) = fld_rhs(i,j,k,bi,bj) -
     &       east      (i,j,k,bi,bj)  * fld_in (i+1, j  , k  ,bi, bj) -
     &       west      (i,j,k,bi,bj)  * fld_out(i-1, j  , k  ,bi, bj) -
     &       north     (i,j,k,bi,bj)  * fld_in (i  , j+1, k  ,bi, bj) -
     &       south     (i,j,k,bi,bj)  * fld_out(i  , j-1, k  ,bi, bj)


            if (k.lt.Nr) fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) -
     &       down      (i,j,k,bi,bj) * fld_in (i  , j  , k+1, bi, bj) -
     &       east_down (i,j,k,bi,bj) * fld_in (i+1, j  , k+1, bi, bj) -
     &       west_down (i,j,k,bi,bj) * fld_in (i-1, j  , k+1, bi, bj) -
     &       north_down(i,j,k,bi,bj) * fld_in (i  , j+1, k+1, bi, bj) -
     &       south_down(i,j,k,bi,bj) * fld_in (i  , j-1, k+1, bi, bj)

            if (k.gt.1) fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) -
     &       up        (i,j,k,bi,bj) * fld_in (i  , j  , k-1 ,bi, bj) -
     &       east_up   (i,j,k,bi,bj) * fld_in (i+1, j  , k-1, bi, bj) -
     &       west_up   (i,j,k,bi,bj) * fld_out(i-1, j  , k-1, bi, bj) -
     &       north_up  (i,j,k,bi,bj) * fld_in (i  , j+1, k-1, bi, bj) -
     &       south_up  (i,j,k,bi,bj) * fld_out(i  , j-1, k-1, bi, bj)

            fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) *
     &          sor*recip_center(i,j,k,bi,bj) +
     &          (oneRL-sor)*fld_in(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo

       ! Update halo region
       CALL EXCH_XYZ_RL ( fld_out , myThid )

       do bj = myByLo(myThid),myByHi(myThid)
        do bi = myBxLo(myThid),myBxHi(myThid)
         do k = 1,Nr
          do j = 1,sNy
           do i = 1,sNx

            ! Compute error: b - Ax
            Ax(i,j,k,bi,bj) =
     &          center(i,j,k,bi,bj) * fld_out(i  , j  , k  , bi, bj) +
     &          east  (i,j,k,bi,bj) * fld_out(i+1, j,   k  , bi, bj) +
     &          west  (i,j,k,bi,bj) * fld_out(i-1, j,   k  , bi, bj) +
     &          north (i,j,k,bi,bj) * fld_out(i  , j+1, k  , bi, bj) +
     &          south (i,j,k,bi,bj) * fld_out(i  , j-1, k  , bi, bj)

            if (k.lt.Nr) Ax(i,j,k,bi,bj) = Ax(i,j,k,bi,bj) +
     &       down      (i,j,k,bi,bj) * fld_out(i  , j  , k+1, bi, bj) +
     &       east_down (i,j,k,bi,bj) * fld_out(i+1, j  , k+1, bi, bj) +
     &       west_down (i,j,k,bi,bj) * fld_out(i-1, j  , k+1, bi, bj) +
     &       north_down(i,j,k,bi,bj) * fld_out(i  , j+1, k+1, bi, bj) +
     &       south_down(i,j,k,bi,bj) * fld_out(i  , j-1, k+1, bi, bj)

            if (k.gt.1 ) Ax(i,j,k,bi,bj) = Ax(i,j,k,bi,bj) +
     &       up        (i,j,k,bi,bj) * fld_out(i  , j  , k-1, bi, bj) +
     &       east_up   (i,j,k,bi,bj) * fld_out(i+1, j  , k-1, bi, bj) +
     &       west_up   (i,j,k,bi,bj) * fld_out(i-1, j  , k-1, bi, bj) +
     &       north_up  (i,j,k,bi,bj) * fld_out(i  , j+1, k-1, bi, bj) +
     &       south_up  (i,j,k,bi,bj) * fld_out(i  , j-1, k-1, bi, bj)

            tile_absErr(bi,bj) = tile_absErr(bi,bj) +
     &          (fld_rhs(i,j,k,bi,bj) - Ax(i,j,k,bi,bj))**2

            ! Reset initial estimate
            fld_in(i,j,k,bi,bj)=fld_out(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo

       CALL GLOBAL_SUM_TILE_RL(tile_absErr, absErr, myThid)

       relErr = absErr / fld_in_norm

       ! Update halos
       CALL EXCH_XYZ_RL ( fld_in , myThid )

      enddo ! end dowhile <TOL and <n_iter_max

      if ( debugLevel.GE.debLevC ) then

        write(msgBuf,'(2A,I10,A,E22.14)') 'SMOOTH_JACOBI_3D REDI: ',
     &      'final iter: ', jacobi_iter,
     &      '  final error: ', SQRT(relErr)
        call print_message(msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid)

        write(msgBuf,'(2A,I10,A,E22.14)') 'SMOOTH_JACOBI_3D REDI: ',
     &      'n_iter_max: ', smooth3DJacobiMaxIters(smoothOpNb),
     &      '  TOL: ', TOL
        call print_message(msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid)
      endif

      RETURN
      END
