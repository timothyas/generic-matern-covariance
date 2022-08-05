#include "SMOOTH_OPTIONS.h"

      SUBROUTINE SMOOTH_READPARMS( myThid )

C     *==========================================================*
C     | SUBROUTINE smooth_readparms
C     | o Routine that reads the pkg/smooth namelist from data.smooth
C     *==========================================================*

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "PARAMS.h"

#include "SMOOTH.h"

c     == routine arguments ==

      integer myThid

C     !FUNCTIONS
      INTEGER ILNBLNK
      EXTERNAL ILNBLNK

c     == local variables ==
C     msgBuf      - Informational/error message buffer
C     iUnit       - Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) namBuf
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit, num_file, num_var, iLen
      integer smoothOpNb

c     == end of interface ==

c--   Read the namelist input.
      namelist /smooth_nml/
     &                   smoothMdsDir,
     &                   smoothMdsDirCreate,
     &                   smooth2Dnbt,
     &                   smooth2DCreateOperator,
     &                   smooth2DCalcNormFactor,
     &                   smooth2DWriteSamples,
     &                   smooth2DConstHorizontal,
     &                   smooth2D_Lx0,
     &                   smooth2D_Ly0,
     &                   smooth2DNbRand,
     &                   smooth2DMaskName,
     &                   smooth2DAlgorithm,
     &                   smooth2DDims,
     &                   smooth2DJacobiMaxIters,
     &                   smooth2DSOROmega,
     &                   smooth2DEigenSize,
C
     &                   smooth3Dnbt,
     &                   smooth3DCreateOperator,
     &                   smooth3DCalcNormFactor,
     &                   smooth3DWriteSamples,
     &                   smooth3DConstHorizontal,
     &                   smooth3DConstVertical,
     &                   smooth3DUseDRCMaxLz,
     &                   smooth3D_Lx0,
     &                   smooth3D_Ly0,
     &                   smooth3D_Lz0,
     &                   smooth3DNbRand,
     &                   smooth3DMaskName,
     &                   smooth3DAlgorithm,
     &                   smooth3DJacobiMaxIters,
     &                   smooth3DSOROmega,
     &                   smooth3DEllipticTol,
     &                   smooth3DNumApplications,
     &                   smooth3DMode


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( .NOT.useSMOOTH ) THEN
C-    pkg SMOOTH is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.smooth is found
         CALL PACKAGES_UNUSED_MSG( 'useSMOOTH', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER( myThid )

c--   Set default values.
      DO smoothOpNb=1,smoothOpNbMax
        smooth2Dnbt(smoothOpNb)             = 0
        smooth2DCreateOperator(smoothOpNb)  = .FALSE.
        smooth2DCalcNormFactor(smoothOpNb)  = .FALSE.
        smooth2DWriteSamples(smoothOpNb)    = .FALSE.
        smooth2DConstHorizontal(smoothOpNb) = .TRUE.
        smooth2D_Lx0(smoothOpNb)            = 0. _d 0
        smooth2D_Ly0(smoothOpNb)            = 0. _d 0
        smooth2DMaskName(smoothOpNb)(1:5)   = 'maskC'
        smooth2DAlgorithm(smoothOpNb)(1:11) = 'WC01       '
        smooth2DDims(smoothOpNb)(1:3)       = 'xy '
        smooth2DNbRand(smoothOpNb)          = 1000
        smooth2DJacobiMaxIters(smoothOpNb)  = 1000
        smooth2DSOROmega(smoothOpNb)        = 1. _d 0
        smooth2DEigenSize(smoothOpNb)       = maxEigenSize

        smooth3Dnbt(smoothOpNb)             = 0
        smooth3DCreateOperator(smoothOpNb)  = .FALSE.
        smooth3DCalcNormFactor(smoothOpNb)  = .FALSE.
        smooth3DWriteSamples(smoothOpNb)    = .FALSE.
        smooth3DConstHorizontal(smoothOpNb) = .TRUE.
        smooth3DConstVertical(smoothOpNb)   = .TRUE.
        smooth3DUseDRCMaxLz(smoothOpNb)     = .FALSE.
        smooth3D_Lx0(smoothOpNb)            = 0. _d 0
        smooth3D_Ly0(smoothOpNb)            = 0. _d 0
        smooth3D_Lz0(smoothOpNb)            = 0. _d 0
        smooth3DMaskName(smoothOpNb)(1:5)   = 'maskC'
        smooth3DAlgorithm(smoothOpNb)(1:11) = 'WC01       '
        smooth3DNbRand(smoothOpNb)          = 1000
        smooth3DJacobiMaxIters(smoothOpNb)  = 1000
        smooth3DSOROmega(smoothOpNb)        = 1. _d 0
        smooth3DEllipticTol(smoothOpNb)     = 1. _d -15
        smooth3DNumApplications(smoothOpNb) = 1
        smooth3DMode(smoothOpNb)            = 'grid'
      ENDDO

      smoothMdsDir       = ' '
      smoothMdsDirCreate = .TRUE.

C--   Read settings from model parameter file "data.smooth".
      WRITE(msgBuf,'(A)') 'SMOOTH_READPARMS: opening data.smooth'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

      CALL OPEN_COPY_DATA_FILE(
     I                          'data.smooth', 'SMOOTH_READPARMS',
     O                          iUnit,
     I                          myThid )

      READ(  iUnit, nml = smooth_nml )

      WRITE(msgBuf,'(2A)') 'SMOOTH_READPARMS: ',
     &       'finished reading data.smooth'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                  SQUEEZE_RIGHT , 1)

#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     remove trailing "/":
      iLen = ILNBLNK( smoothMdsDir )
      IF ( iLen.GE.2 ) THEN
       IF ( smoothMdsDir(iLen:iLen).EQ.'/' ) THEN
         namBuf = smoothMdsDir
         WRITE(smoothMdsDir,'(A)') namBuf(1:iLen-1)
       ENDIF
      ENDIF

C--   Print pkg/smooth settings to standard output:
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      WRITE(msgBuf,'(A)') '// pkg/smooth configuration'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &  SQUEEZE_RIGHT , 1)

      DO smoothOpNb=1,smoothOpNbMax
       WRITE(msgBuf,'(A,I2,I6,2f6.0,4A,A,L)') 'smooth 2D parameters: ',
     &    smoothOpNb,smooth2Dnbt(smoothOpNb),
     &    smooth2D_Lx0(smoothOpNb),smooth2D_Ly0(smoothOpNb),
     &    '  ', smooth2DmaskName(smoothOpNb),
     &    '  ', smooth2DAlgorithm(smoothOpNb),
     &    '  ', smooth2DCalcNormFactor(smoothOpNb)
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT , 1)
      ENDDO

      DO smoothOpNb=1,smoothOpNbMax
       WRITE(msgBuf,'(A,I2,I6,3f6.0,4A)') 'smooth 3D parameters: ',
     &    smoothOpNb,smooth3Dnbt(smoothOpNb),
     &    smooth3D_Lx0(smoothOpNb),smooth3D_Ly0(smoothOpNb),
     &    smooth3D_Lz0(smoothOpNb),
     &    '  ', smooth3DmaskName(smoothOpNb),
     &    '  ', smooth3DAlgorithm(smoothOpNb)
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT , 1)
      ENDDO

      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      WRITE(msgBuf,'(A)') '// End of pkg/smooth config. summary'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &  SQUEEZE_RIGHT , 1)
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &  SQUEEZE_RIGHT , 1)


      _END_MASTER( myThid )

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

      RETURN
      END
