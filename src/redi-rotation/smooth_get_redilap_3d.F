#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      SUBROUTINE SMOOTH_GET_REDILAP_3D (
     U     center, east, west, north, south, up, down,
     U     east_up, east_down, west_up, west_down,
     U     north_up, north_down, south_up, south_down,
     I     smoothOpNb, smooth3DMask, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_get_redilap_3d
C     | o get the matrix to the discretized 3D
C     |   laplacian operator: div( kappa * grad( ) )
C     | o Use Redi elements of tensor:
C     |
C     |             Kux Kuy Kuz    Kux  0  Kuz
C     |    kappa =  Kvx Kvy Kvz ->  0  Kvy Kvz
C     |             Kwx Kwy Kwz    Kwx Kwy Kwz
C     |
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
      _RL center    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL east_up   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east_down (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west_up   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west_down (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north_up  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north_down(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south_up  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south_down(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      INTEGER smoothOpNb
      _RS smooth3Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER myThid

C --- Local parameters
      _RL xA (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL yA (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL zAx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL zAy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL m1mask, p1mask
      _RL denom_west, denom_east
      _RL denom_south, denom_north
      _RL denom_wxu, denom_wyu
      _RL denom_wxd, denom_wyd
      INTEGER i,j,k,bi,bj
      INTEGER km1, kp1
      CHARACTER*(MAX_LEN_FNAM) fnamegeneric

C --- Get diagonal elements
      CALL SMOOTH_GET_LAPLACIAN_3D(
     &      center, east, west, north, south, up, down,
     &      smoothOpNb, smooth3DMask, myThid)

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DE1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, east, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DW1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, west, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DN1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, north, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DS1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, south, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DC1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, center, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DU1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, up, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DD1',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, down, 1, 1, myThid )

C --- Initialize output
      DO bj = myByLo(myThid),myByHi(myThid)
       DO bi = myBxLo(myThid),myBxHi(myThid)
        DO k = 1,Nr
         km1 = MAX(k-1, 1)
         kp1 = MIN(k+1, Nr)
         m1mask = oneRL
         p1mask = oneRL
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx
           east_up   (i,j,k,bi,bj) = zeroRL
           east_down (i,j,k,bi,bj) = zeroRL
           west_up   (i,j,k,bi,bj) = zeroRL
           west_down (i,j,k,bi,bj) = zeroRL
           north_up  (i,j,k,bi,bj) = zeroRL
           north_down(i,j,k,bi,bj) = zeroRL
           south_up  (i,j,k,bi,bj) = zeroRL
           south_down(i,j,k,bi,bj) = zeroRL

           xA(i,j,k,bi,bj) = _dyG(i,j,bi,bj) * drF(k) *
     &         _hFacW(i,j,k,bi,bj) * smooth3D_Kuz(i,j,k,bi,bj) *
     &         smooth3DMask(i-1, j  , k  , bi, bj) *
     &         smooth3DMask(i  , j  , k  , bi, bj)

           yA(i,j,k,bi,bj) = _dxG(i,j,bi,bj) * drF(k) *
     &         _hFacS(i,j,k,bi,bj) * smooth3D_Kvz(i,j,k,bi,bj) *
     &         smooth3DMask(i  , j-1, k  , bi, bj) *
     &         smooth3DMask(i  , j  , k  , bi, bj)

           zAx(i,j,k,bi,bj) = _rA(i,j,bi,bj) *
     &         _hFacC(i,j,k,bi,bj) * smooth3D_Kwx(i,j,k,bi,bj) *
     &         smooth3DMask(i  , j  , km1, bi, bj) *
     &         smooth3DMask(i  , j  , k  , bi, bj) *
     &         m1mask

           zAy(i,j,k,bi,bj) = _rA(i,j,bi,bj) *
     &         _hFacC(i,j,k,bi,bj) * smooth3D_Kwy(i,j,k,bi,bj) *
     &         smooth3DMask(i  , j  , km1, bi, bj) *
     &         smooth3DMask(i  , j  , k  , bi, bj) *
     &         m1mask
          ENDDO
         ENDDO
        ENDDO

        DO k = 1,Nr

C-       These are purely so I don't have to write if statements for
C        k==1 or k==Nr
         km1 = MAX(k-1, 1)
         kp1 = MIN(k+1, Nr)
         m1mask = oneRL
         p1mask = oneRL
         IF (k.EQ.1 ) m1mask = zeroRL
         IF (k.EQ.Nr) p1mask = zeroRL

         DO j = 1-OLy+1,sNy+OLy-1
          DO i = 1-OLx+1,sNx+OLx-1


C-        Useful
          denom_west  = oneRL /
     &      (rA (i, j, bi, bj) + rA (i-1, j  , bi, bj)) /
     &      (drC(k) + drC(kp1))
          denom_east  = oneRL /
     &      (rA (i, j, bi, bj) + rA (i+1, j  , bi, bj)) /
     &      (drC(k) + drC(kp1))
          denom_south = oneRL /
     &      (rA (i, j, bi, bj) + rA (i  , j-1, bi, bj)) /
     &      (drC(k) + drC(kp1))
          denom_north = oneRL /
     &      (rA (i, j, bi, bj) + rA (i  , j+1, bi, bj)) /
     &      (drC(k) + drC(kp1))

          denom_wxu    = oneRL /
     &      (dyG(i, j, bi, bj) + dyG(i+1, j  , bi, bj)) /
     &      (drF(k) + drF(km1))

          denom_wxd    = oneRL /
     &      (dyG(i, j, bi, bj) + dyG(i+1, j  , bi, bj)) /
     &      (drF(k) + drF(kp1))

          denom_wyu    = oneRL /
     &      (dxG(i, j, bi, bj) + dxG(i  , j+1, bi, bj)) /
     &      (drF(k) + drF(km1))

          denom_wyd    = oneRL /
     &      (dxG(i, j, bi, bj) + dxG(i  , j+1, bi, bj)) /
     &      (drF(k) + drF(kp1))


C-        West
          west_up(i,j,k,bi,bj) = west_up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_west * m1mask *
     &      rA(i-1, j, bi, bj) * drC(k) *
     &      smooth3Dmask(i-1,j,k,bi,bj)*smooth3Dmask(i-1,j,km1,bi,bj) *
     &      xA(i,j,k,bi,bj)

          west_down(i,j,k,bi,bj) = west_down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_west * p1mask *
     &      rA(i-1, j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i-1,j,k,bi,bj)*smooth3Dmask(i-1,j,kp1,bi,bj) *
     &      xA(i,j,k,bi,bj)

          west(i,j,k,bi,bj) = west(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_west * m1mask *
     &      rA(i-1, j, bi, bj) * drC(k) *
     &      smooth3Dmask(i-1,j,k,bi,bj)*smooth3Dmask(i-1,j,km1,bi,bj) *
     &      xA(i,j,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_west * p1mask *
     &      rA(i-1, j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i-1,j,k,bi,bj)*smooth3Dmask(i-1,j,kp1,bi,bj) *
     &      xA(i,j,k,bi,bj)

          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_west * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,km1,bi,bj) *
     &      xA(i,j,k,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_west * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,kp1,bi,bj) *
     &      xA(i,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_west * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,km1,bi,bj) *
     &      xA(i,j,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_west * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,kp1,bi,bj) *
     &      xA(i,j,k,bi,bj)

C-        East
          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_east * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,km1,bi,bj) *
     &      xA(i+1,j,k,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_east * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,kp1,bi,bj) *
     &      xA(i+1,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_east * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,km1,bi,bj) *
     &      xA(i+1,j,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_east * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i  ,j,k,bi,bj)*smooth3Dmask(i  ,j,kp1,bi,bj) *
     &      xA(i+1,j,k,bi,bj)

          east_up(i,j,k,bi,bj) = east_up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_east * m1mask *
     &      rA(i+1, j, bi, bj) * drC(k) *
     &      smooth3Dmask(i+1,j,k,bi,bj)*smooth3Dmask(i+1,j,km1,bi,bj) *
     &      xA(i+1,j,k,bi,bj)

          east_down(i,j,k,bi,bj) = east_down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_east * p1mask *
     &      rA(i+1, j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i+1,j,k,bi,bj)*smooth3Dmask(i+1,j,kp1,bi,bj) *
     &      xA(i+1,j,k,bi,bj)

          east(i,j,k,bi,bj) = east(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_east * m1mask *
     &      rA(i+1, j, bi, bj) * drC(k) *
     &      smooth3Dmask(i+1,j,k,bi,bj)*smooth3Dmask(i+1,j,km1,bi,bj) *
     &      xA(i+1,j,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_east * p1mask *
     &      rA(i+1, j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i+1,j,k,bi,bj)*smooth3Dmask(i+1,j,kp1,bi,bj) *
     &      xA(i+1,j,k,bi,bj)


C-        South
          south_up(i,j,k,bi,bj) = south_up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_south * m1mask *
     &      rA(i, j-1, bi, bj) * drC(k) *
     &      smooth3Dmask(i,j-1,k,bi,bj)*smooth3Dmask(i,j-1,km1,bi,bj) *
     &      yA(i,j,k,bi,bj)

          south_down(i,j,k,bi,bj) = south_down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_south * p1mask *
     &      rA(i, j-1, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j-1,k,bi,bj)*smooth3Dmask(i,j-1,kp1,bi,bj) *
     &      yA(i,j,k,bi,bj)

          south(i,j,k,bi,bj) = south(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_south * m1mask *
     &      rA(i, j-1, bi, bj) * drC(k) *
     &      smooth3Dmask(i,j-1,k,bi,bj)*smooth3Dmask(i,j-1,km1,bi,bj) *
     &      yA(i,j,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_south * p1mask *
     &      rA(i, j-1, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j-1,k,bi,bj)*smooth3Dmask(i,j-1,kp1,bi,bj) *
     &      yA(i,j,k,bi,bj)

          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_south * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,km1,bi,bj) *
     &      yA(i,j,k,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_south * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,kp1,bi,bj) *
     &      yA(i,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_south * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,km1,bi,bj) *
     &      yA(i,j,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_south * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,kp1,bi,bj) *
     &      yA(i,j,k,bi,bj)


C-        North
          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_north * m1mask *
     &      rA(i, j  , bi, bj) * drC(k) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,km1,bi,bj) *
     &      yA(i,j+1,k,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_north * p1mask *
     &      rA(i, j  , bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,kp1,bi,bj) *
     &      yA(i,j+1,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_north * m1mask *
     &      rA(i, j  , bi, bj) * drC(k) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,km1,bi,bj) *
     &      yA(i,j+1,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_north * p1mask *
     &      rA(i, j  , bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j  ,k,bi,bj)*smooth3Dmask(i,j  ,kp1,bi,bj) *
     &      yA(i,j+1,k,bi,bj)

          north_up(i,j,k,bi,bj) = north_up(i,j,k,bi,bj) +
     &      recip_drC(k) * denom_north * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i,j+1,k,bi,bj)*smooth3Dmask(i,j+1,km1,bi,bj) *
     &      yA(i,j+1,k,bi,bj)

          north_down(i,j,k,bi,bj) = north_down(i,j,k,bi,bj) +
     &      -recip_drC(kp1) * denom_north * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j+1,k,bi,bj)*smooth3Dmask(i,j+1,kp1,bi,bj) *
     &      yA(i,j+1,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_drC(k) * denom_north * m1mask *
     &      rA(i  , j, bi, bj) * drC(k) *
     &      smooth3Dmask(i,j+1,k,bi,bj)*smooth3Dmask(i,j+1,km1,bi,bj) *
     &      yA(i,j+1,k,bi,bj) +
C
     &      recip_drC(kp1) * denom_north * p1mask *
     &      rA(i  , j, bi, bj) * drC(kp1) *
     &      smooth3Dmask(i,j+1,k,bi,bj)*smooth3Dmask(i,j+1,kp1,bi,bj) *
     &      yA(i,j+1,k,bi,bj)


C-        Top, X
          west_up(i,j,k,bi,bj) = west_up(i,j,k,bi,bj) +
     &      -recip_dxC(i  ,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i  ,j, bi, bj) * drF(km1) *
     &      smooth3Dmask(i-1,j,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      recip_dxC(i  ,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i  ,j, bi, bj) * drF(km1) *
     &      smooth3Dmask(i-1,j,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          west(i,j,k,bi,bj) = west(i,j,k,bi,bj) +
     &      -recip_dxC(i  ,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i  ,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i-1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      recip_dxC(i  ,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i  ,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i-1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          east_up(i,j,k,bi,bj) = east_up(i,j,k,bi,bj) +
     &      recip_dxC(i+1,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i+1,j, bi, bj) * drF(km1) *
     &      smooth3Dmask(i+1,j,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      -recip_dxC(i+1,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i+1,j, bi, bj) * drF(km1) *
     &      smooth3Dmask(i+1,j,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          east(i,j,k,bi,bj) = east(i,j,k,bi,bj) +
     &      recip_dxC(i+1,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i+1,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i+1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_dxC(i+1,j,bi,bj) * denom_wxu * m1mask *
     &      dyG(i+1,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i+1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,k,bi,bj)


C-        Top, Y
          south_up(i,j,k,bi,bj) = south_up(i,j,k,bi,bj) +
     &      -recip_dyC(i,j  ,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j  , bi, bj) * drF(km1) *
     &      smooth3Dmask(i,j-1,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      recip_dyC(i,j,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j  , bi, bj) * drF(km1) *
     &      smooth3Dmask(i,j-1,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          south(i,j,k,bi,bj) = south(i,j,k,bi,bj) +
     &      -recip_dyC(i,j,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j  , bi, bj) * drF(k) *
     &      smooth3Dmask(i,j-1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      recip_dyC(i,j,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j  , bi, bj) * drF(k) *
     &      smooth3Dmask(i,j-1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          north_up(i,j,k,bi,bj) = north_up(i,j,k,bi,bj) +
     &      recip_dyC(i,j+1,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j+1, bi, bj) * drF(km1) *
     &      smooth3Dmask(i,j+1,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          up(i,j,k,bi,bj) = up(i,j,k,bi,bj) +
     &      -recip_dyC(i,j+1,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j+1, bi, bj) * drF(km1) *
     &      smooth3Dmask(i,j+1,km1,bi,bj)*smooth3Dmask(i,j,km1,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          north(i,j,k,bi,bj) = north(i,j,k,bi,bj) +
     &      recip_dyC(i,j+1,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j+1, bi, bj) * drF(k) *
     &      smooth3Dmask(i,j+1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,k,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_dyC(i,j+1,bi,bj) * denom_wyu * m1mask *
     &      dxG(i,j+1, bi, bj) * drF(k) *
     &      smooth3Dmask(i,j+1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,k,bi,bj)

C-        Bottom, X
          west(i,j,k,bi,bj) = west(i,j,k,bi,bj) +
     &      -recip_dxC(i  ,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i  ,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i-1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      recip_dxC(i  ,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i  ,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i-1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          west_down(i,j,k,bi,bj) = west_down(i,j,k,bi,bj) +
     &      -recip_dxC(i  ,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i  ,j, bi, bj) * drF(kp1) *
     &      smooth3Dmask(i-1,j,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      recip_dxC(i  ,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i  ,j, bi, bj) * drF(kp1) *
     &      smooth3Dmask(i-1,j,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          east(i,j,k,bi,bj) = east(i,j,k,bi,bj) +
     &      recip_dxC(i+1,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i+1,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i+1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_dxC(i+1,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i+1,j, bi, bj) * drF(k) *
     &      smooth3Dmask(i+1,j,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          east_down(i,j,k,bi,bj) = east_down(i,j,k,bi,bj) +
     &      recip_dxC(i+1,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i+1,j, bi, bj) * drF(kp1) *
     &      smooth3Dmask(i+1,j,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      -recip_dxC(i+1,j,bi,bj) * denom_wxd * p1mask *
     &      dyG(i+1,j, bi, bj) * drF(kp1) *
     &      smooth3Dmask(i+1,j,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAx(i,j,kp1,bi,bj)


C-        Bottom, Y
          south(i,j,k,bi,bj) = south(i,j,k,bi,bj) +
     &      -recip_dyC(i,j  ,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j  , bi, bj) * drF(k) *
     &      smooth3Dmask(i,j-1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      recip_dyC(i,j,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j  , bi, bj) * drF(k) *
     &      smooth3Dmask(i,j-1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          south_down(i,j,k,bi,bj) = south_down(i,j,k,bi,bj) +
     &      -recip_dyC(i,j,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j  , bi, bj) * drF(kp1) *
     &      smooth3Dmask(i,j-1,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      recip_dyC(i,j,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j  , bi, bj) * drF(kp1) *
     &      smooth3Dmask(i,j-1,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          north(i,j,k,bi,bj) = north(i,j,k,bi,bj) +
     &      recip_dyC(i,j+1,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j+1, bi, bj) * drF(k) *
     &      smooth3Dmask(i,j+1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          center(i,j,k,bi,bj) = center(i,j,k,bi,bj) +
     &      -recip_dyC(i,j+1,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j+1, bi, bj) * drF(k) *
     &      smooth3Dmask(i,j+1,k  ,bi,bj)*smooth3Dmask(i,j,k  ,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          north_down(i,j,k,bi,bj) = north_down(i,j,k,bi,bj) +
     &      recip_dyC(i,j+1,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j+1, bi, bj) * drF(kp1) *
     &      smooth3Dmask(i,j+1,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          down(i,j,k,bi,bj) = down(i,j,k,bi,bj) +
     &      -recip_dyC(i,j+1,bi,bj) * denom_wyd * p1mask *
     &      dxG(i,j+1, bi, bj) * drF(kp1) *
     &      smooth3Dmask(i,j+1,kp1,bi,bj)*smooth3Dmask(i,j,kp1,bi,bj) *
     &      zAy(i,j,kp1,bi,bj)

          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XYZ_RL( center, myThid )
      CALL EXCH_XYZ_RL( east, myThid )
      CALL EXCH_XYZ_RL( west, myThid )
      CALL EXCH_XYZ_RL( north, myThid )
      CALL EXCH_XYZ_RL( south, myThid )
      CALL EXCH_XYZ_RL( up, myThid )
      CALL EXCH_XYZ_RL( down, myThid )

      CALL EXCH_XYZ_RL( east_up, myThid )
      CALL EXCH_XYZ_RL( east_down, myThid )
      CALL EXCH_XYZ_RL( west_up, myThid )
      CALL EXCH_XYZ_RL( west_down, myThid )
      CALL EXCH_XYZ_RL( north_up, myThid )
      CALL EXCH_XYZ_RL( north_down, myThid )
      CALL EXCH_XYZ_RL( south_up, myThid )
      CALL EXCH_XYZ_RL( south_down, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DEUp',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, east_up, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DEDn',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, east_down, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DWUp',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, west_up, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DWDn',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, west_down, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DNUp',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, north_up, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DNDn',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, north_down, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DSUp',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, south_up, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DSDn',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, south_down, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DE',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, east, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DW',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, west, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DN',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, north, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DS',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, south, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DC',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, center, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DU',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, up, 1, 1, myThid )

      WRITE(fnamegeneric(1:MAX_LEN_FNAM),'(1a,i3.3)')
     &   'smooth3DD',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, down, 1, 1, myThid )
      RETURN
      END
