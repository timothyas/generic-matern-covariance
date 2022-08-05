#include "SMOOTH_OPTIONS.h"

      SUBROUTINE SMOOTH_INVERSE_MATERN_3D (
     U     fld_in,
     I     smooth3DMask,smoothOpNb,doNormalize,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_inverse_matern_3d
C     | o Given the laplace-like differential operator
C     |
C     |         A=( delta - div( kappa * grad() )
C     |
C     |   this employs an iterative method to find the solution
C     |
C     |         u = A^{-1}b
C     |
C     | o This is used as a preconditioner to get an estimate
C     |   of the prior preconditioned misfit Hessian (Gauss Newton
C     |   approximation)
C     |
C     | o if doNormalize==True, result is multiplied by 1/sqrt(var(A^{-1}))
C     |   which is written to file as smooth2Dnorm00<smoothOpNb>
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      _RL fld_in      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER smoothOpNb
      LOGICAL doNormalize
      INTEGER myThid

C --- Local parameters
      INTEGER i,j,k,bi,bj
      INTEGER app
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL east_up   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west_up   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north_up  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south_up  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east_down (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west_down (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north_down(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south_down(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

C ------------------------------------------------------------------
C the actual algorithm
C ------------------------------------------------------------------

c --- 0. Apply inv sqrt of cell volume and determinant of mapping:
      if (doNormalize) then
       do bj = myByLo(myThid),myByHi(myThid)
        do bi = myBxLo(myThid),myBxHi(myThid)
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            fld_in(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)*
     &          smooth3DRandNorm(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo
      endif

c --- 1. Get laplacian matrix elements
      IF (smooth3DMode(smoothOpNb).EQ.'grid') THEN
       call smooth_get_laplacian_3d(
     &       center, east, west, north, south, up, down,
     &       smoothOpNb, smooth3DMask, myThid )
      ELSEIF (smooth3DMode(smoothOpNb).EQ.'redi') THEN
       call smooth_get_redilap_3d(
     &       center, east, west, north, south, up, down,
     &       east_up, east_down, west_up, west_down,
     &       north_up, north_down, south_up, south_down,
     &       smoothOpNb, smooth3DMask, myThid )
      ENDIF


c --- 2. Create A matrix associated with operator A = (delta - nabla^2)
c        operate directly on overlap region, inputs "pre-exchanged"
      do bj = myByLo(myThid),myByHi(myThid)
       do bi = myBxLo(myThid),myBxHi(myThid)
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           center(i,j,k,bi,bj) = -center(i,j,k,bi,bj) +
     &          smooth3DDelta(i,j,k,bi,bj)
           east (i,j,k,bi,bj)  = -east (i,j,k,bi,bj)
           west (i,j,k,bi,bj)  = -west (i,j,k,bi,bj)
           north(i,j,k,bi,bj)  = -north(i,j,k,bi,bj)
           south(i,j,k,bi,bj)  = -south(i,j,k,bi,bj)
           up   (i,j,k,bi,bj)  = -up   (i,j,k,bi,bj)
           down (i,j,k,bi,bj)  = -down (i,j,k,bi,bj)
           IF (smooth3DMode(smoothOpNb).EQ.'redi') THEN
            east_up   (i,j,k,bi,bj) = -east_up   (i,j,k,bi,bj)
            west_up   (i,j,k,bi,bj) = -west_up   (i,j,k,bi,bj)
            north_up  (i,j,k,bi,bj) = -north_up  (i,j,k,bi,bj)
            south_up  (i,j,k,bi,bj) = -south_up  (i,j,k,bi,bj)
            east_down (i,j,k,bi,bj) = -east_down (i,j,k,bi,bj)
            west_down (i,j,k,bi,bj) = -west_down (i,j,k,bi,bj)
            north_down(i,j,k,bi,bj) = -north_down(i,j,k,bi,bj)
            south_down(i,j,k,bi,bj) = -south_down(i,j,k,bi,bj)
           ENDIF
          enddo
         enddo
        enddo
       enddo
      enddo

c --- 3. Get jacobi to solve our problem
      do app = 1, smooth3DNumApplications(smoothOpNb)
       IF (smooth3DMode(smoothOpNb).EQ.'grid') THEN
        call smooth_jacobi_3d( fld_in, smoothOpNb,
     &        center, east, west, north, south, up, down,
     &        smooth3DMask, myThid )
       ELSEIF (smooth3DMode(smoothOpNb).EQ.'redi') THEN
        call smooth_sor_redi_3d( fld_in, smoothOpNb,
     &        center, east, west, north, south, up, down,
     &        east_up, east_down, west_up, west_down,
     &        north_up, north_down, south_up, south_down,
     &        smooth3DMask, myThid )
       ENDIF
      enddo

c --- 4. Apply normalization field [i.e. 1/sqrt(var(filter))]:
      if (doNormalize) then
       do bj = myByLo(myThid),myByHi(myThid)
        do bi = myBxLo(myThid),myBxHi(myThid)
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            fld_in(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)*
     &          smooth3DNorm(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo
      endif

      RETURN
      END
