C pkg/smooth constants

      INTEGER     smoothprec
      PARAMETER ( smoothprec = 32 )

      LOGICAL smooth3DdoImpldiff
      PARAMETER ( smooth3DdoImpldiff = .TRUE. )

      INTEGER smoothOpNbMax
      PARAMETER ( smoothOpNbMax = 10 )

      _RL smooth2DdelTime, smooth3DdelTime
      PARAMETER ( smooth2DdelTime = 1. _d 0 )
      PARAMETER ( smooth3DdelTime = 1. _d 0 )

      INTEGER maxEigenSize
      PARAMETER ( maxEigenSize = 400 )

C parameters:

      COMMON /smooth_operators_i/
     & smooth2Dnbt,
     & smooth2DNbRand,
     & smooth2DJacobiMaxIters,
     & smooth2DEigenSize,
     & smooth3Dnbt,
     & smooth3DNbRand,
     & smooth3DJacobiMaxIters,
     & smooth3DNumApplications
      INTEGER smooth2Dnbt(smoothOpNbMax)
      INTEGER smooth2DNbRand(smoothOpNbMax)
      INTEGER smooth2DJacobiMaxIters(smoothOpNbMax)
      INTEGER smooth2DEigenSize(smoothOpNbMax)
      INTEGER smooth3Dnbt(smoothOpNbMax)
      INTEGER smooth3DNbRand(smoothOpNbMax)
      INTEGER smooth3DJacobiMaxIters(smoothOpNbMax)
      INTEGER smooth3DNumApplications(smoothOpNbMax)

      COMMON /smooth_param_r/
     & smooth3DtotTime,
     & smooth3D_Lx0, smooth3D_Ly0, smooth3D_Lz0,
     & smooth2DtotTime, smooth2D_Lx0, smooth2D_Ly0,
     & smooth3DSOROmega, smooth2DSOROmega,
     & smooth3DEllipticTol
      _RL smooth3DtotTime,
     & smooth3D_Lx0(smoothOpNbMax),
     & smooth3D_Ly0(smoothOpNbMax), smooth3D_Lz0(smoothOpNbMax)
      _RL smooth2DtotTime,
     & smooth2D_Lx0(smoothOpNbMax), smooth2D_Ly0(smoothOpNbMax)
      _RL smooth3DSOROmega(smoothOpNbMax)
      _RL smooth2DSOROmega(smoothOpNbMax)
      _RL smooth3DEllipticTol(smoothOpNbMax)

      COMMON /smooth_flds_c/
     & smoothMdsDir,
     & smooth3DmaskName, smooth2DmaskName,
     & smooth3DAlgorithm,
     & smooth2DAlgorithm,
     & smooth2DDims,
     & smooth3DMode
      CHARACTER*(MAX_LEN_FNAM) smoothMdsDir
      CHARACTER*(5)  smooth3DmaskName(smoothOpNbMax)
      CHARACTER*(5)  smooth2DmaskName(smoothOpNbMax)
      CHARACTER*(11) smooth3DAlgorithm(smoothOpNbMax)
      CHARACTER*(11) smooth2DAlgorithm(smoothOpNbMax)
      CHARACTER*(3)  smooth2DDims(smoothOpNbMax)
      CHARACTER*(4)  smooth3DMode(smoothOpNbMax)

      COMMON /smooth_param_l/
     & smoothMdsDirCreate,
     & smooth2DConstHorizontal,
     & smooth2DCreateOperator,
     & smooth2DCalcNormFactor,
     & smooth2DWriteSamples,
     & smooth3DConstHorizontal,
     & smooth3DConstVertical,
     & smooth3DCreateOperator,
     & smooth3DCalcNormFactor,
     & smooth3DWriteSamples,
     & smooth3DUseDRCMaxLz
      LOGICAL smoothMdsDirCreate
      LOGICAL smooth2DConstHorizontal(smoothOpNbMax)
      LOGICAL smooth2DCreateOperator (smoothOpNbMax)
      LOGICAL smooth2DCalcNormFactor (smoothOpNbMax)
      LOGICAL smooth2DWriteSamples   (smoothOpNbMax)
      LOGICAL smooth3DConstHorizontal(smoothOpNbMax)
      LOGICAL smooth3DConstVertical  (smoothOpNbMax)
      LOGICAL smooth3DCreateOperator (smoothOpNbMax)
      LOGICAL smooth3DCalcNormFactor (smoothOpNbMax)
      LOGICAL smooth3DWriteSamples   (smoothOpNbMax)
      LOGICAL smooth3DUseDRCMaxLz    (smoothOpNbMax)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C fields:
      COMMON /smooth_flds_rs/
     & smooth_recip_hFacC, smooth_hFacW, smooth_hFacS
      _RS
     & smooth_recip_hFacC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth_hFacW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth_hFacS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      COMMON /smooth_flds_rl/
     & smooth3D_Lx, smooth3D_Ly, smooth3D_Lz,
     & smooth3Dnorm,
     & smooth2D_Lx, smooth2D_Ly,
     & smooth2Dnorm,
     & smoothXZNorm,
     & smoothYZNorm
      _RL
     & smooth3D_Lx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Ly(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Lz(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3Dnorm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL
     & smooth2D_Lx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smooth2D_Ly(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smooth2Dnorm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL
     & smoothXZNorm (1-OLx:sNx+OLx,Nr,nSx,nSy),
     & smoothYZNorm (1-OLy:sNy+OLy,Nr,nSx,nSy)

      COMMON /smooth_operators_r/
     & smooth3D_kappaR, smooth3D_Kwz, smooth3D_Kux, smooth3D_Kvy,
     & smooth3D_Kwx, smooth3D_Kwy, smooth3D_Kuy, smooth3D_Kuz,
     & smooth3D_Kvx, smooth3D_Kvz,
     & smooth3DDelta, smooth3DRandNorm,
     & smooth2D_Kux, smooth2D_Kvy,
     & smooth2DDelta, smooth2DRandNorm,
     & smoothXZ_Kux, smoothXZ_Kwz,
     & smoothXZDelta, smoothXZRandNorm,
     & smoothYZ_Kvy, smoothYZ_Kwz,
     & smoothYZDelta, smoothYZRandNorm
      _RL
     & smooth3D_kappaR (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kwz    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kwx    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kwy    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kux    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kuy    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kuz    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kvy    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kvx    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kvz    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3DDelta   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3DRandNorm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL
     & smooth2D_Kux    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smooth2D_Kvy    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smooth2DDelta   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smooth2DRandNorm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL
     & smoothXZ_Kux    (1-Olx:sNx+Olx,Nr,nSx,nSy),
     & smoothXZ_Kwz    (1-Olx:sNx+Olx,Nr,nSx,nSy),
     & smoothXZDelta   (1-OLx:sNx+OLx,Nr,nSx,nSy),
     & smoothXZRandNorm(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL
     & smoothYZ_Kvy    (1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smoothYZ_Kwz    (1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smoothYZDelta   (1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smoothYZRandNorm(1-OLy:sNy+OLy,Nr,nSx,nSy)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
