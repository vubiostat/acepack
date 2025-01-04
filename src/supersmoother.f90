  !----------------------------------------------------------------------------
 !
! This file is part of acepack.
!
! Copyright 1985,2007 Jerome H. Friedman
! Copyright 2016,2025 Shawn Garbett, Vanderbilt University Medical Center
!
! Permission to use, copy, modify, distribute, and sell this software and
! its documentation for any purpose is hereby granted without fee,
! provided that the above copyright notice appear in all copies and that
! both that copyright notice and this permission notice appear in
! supporting documentation. No representations are made about the
! suitability of this software for any purpose.  It is provided "as is"
! without express or implied warranty.
!______________________________________________________________________________

! Original: https://stacks.stanford.edu/file/druid:gw754yg8889/ORION%20003.pdf
! J Friedman, W Stuetzle. Smoothing of Scatterplots, Stanford Project Orion, July 1982
SUBROUTINE SuperSmoother(X,Y,W,SPAN,DOF,N,CROSS,SMO,S0,RSS,SCRAT)
  IMPLICIT NONE
  
      integer I,IBNEW,IBOLD,IS2,ITNEW,ITOLD,J,JJ,M0,NTIE,N
      double precision  X(N),Y(N),W(N),SMO(N),SCRAT(N),RSS,SPAN
      double precision  XIN, YIN, WIN, XOUT, YOUT
      double precision  WT, ISPAN, DOF
      DOUBLE PRECISION SUMW,XBAR,YBAR,COV,VAR, S0, R
      INTEGER FIXEDS,CROSS,LINE
      LINE=1
      FIXEDS=1
      IF(SPAN .GE. 1.0)GOTO 10131
      LINE=0
10131 CONTINUE
      XBAR=X(1)
      YBAR=Y(1)
      COV=0.
      VAR=0.
      SUMW=W(1)
      IF(LINE .NE. 1)GOTO 10151
      DO 10161 I=2,N
      XIN=X(I)
      YIN=Y(I)
      WIN=W(I)
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW
      SUMW=SUMW+WIN
10161 CONTINUE
      CONTINUE
      I=1
      GOTO 10173
10171 I=I+1
10173 IF((I).GT.(N))GOTO 10172
      IF(.NOT.(CROSS.eq.1))GOTO 10191
      XOUT=X(I)
      YOUT=Y(I)
      WIN=W(I)
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)
      SUMW=SUMW-WIN
10191 CONTINUE
      IF(VAR .LE. 0.)GOTO 10211
      SMO(I)=COV*(X(I)-XBAR)/VAR
      GOTO 10221
10211 CONTINUE
      SMO(I)=0
10221 CONTINUE
      CONTINUE
      IF(.NOT.(CROSS.eq.1))GOTO 10241
      XIN=X(I)
      YIN=Y(I)
      WIN=W(I)
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW
      SUMW=SUMW+WIN
10241 CONTINUE
      GOTO 10171
10172 CONTINUE
      S0=YBAR
      SCRAT(1)=COV/VAR
      DOF=1.0
      GOTO 10251
10151 CONTINUE
      ITOLD=1
      IBOLD=1
      DOF=-1.0
      DO 10261 I=1,N
      SCRAT(I)=Y(I)
10261 CONTINUE
      CONTINUE
      IF(.NOT.(cross.eq.0))GOTO 10281
      I=0
10291 IF(I.GE.N-1) GOTO 10292
      I=I+1
      M0=I
10301 IF(X(I+1).GT.X(I)) GOTO 10302
      I=I+1
      IF(I .LT. N)GOTO 10301
10302 CONTINUE
      IF(I.EQ.M0)GOTO 10291
      NTIE=I-M0+1
      R=0.
      WT=0.
      DO 10311 JJ=M0,I
      J=JJ
      R=R+Y(J)*W(J)
      WT=WT+W(J)
10311 CONTINUE
      CONTINUE
      R=R/WT
      DO 10321 J=M0,I
      Y(J)=R
10321 CONTINUE
      CONTINUE
      GOTO 10291
10292 CONTINUE
10281 CONTINUE
      ISPAN=N*SPAN
      IF(.NOT.(FIXEDS.eq.1))GOTO 10341
      IS2=INT(ISPAN/2)
      IF(IS2 .GE. 1)GOTO 10361
      IS2=1
10361 CONTINUE
10341 CONTINUE
      DO 10371 I=1,N
      ITNEW=MIN(I+IS2,N)
      IBNEW=MAX(I-IS2,1)
10381 IF(ITOLD .GE. ITNEW) GOTO 10382
      ITOLD=ITOLD+1
      XIN=X(ITOLD)
      YIN=Y(ITOLD)
      WIN=W(ITOLD)
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW
      SUMW=SUMW+WIN
      GOTO 10381
10382 CONTINUE
10391 IF(IBOLD .LE. IBNEW) GOTO 10392
      IBOLD=IBOLD-1
      XIN=X(IBOLD)
      YIN=Y(IBOLD)
      WIN=W(IBOLD)
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW
      SUMW=SUMW+WIN
      GOTO 10391
10392 CONTINUE
10401 IF(ITOLD .LE. ITNEW) GOTO 10402
      XOUT=X(ITOLD)
      YOUT=Y(ITOLD)
      WIN=W(ITOLD)
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)
      SUMW=SUMW-WIN
      ITOLD=ITOLD-1
      GOTO 10401
10402 CONTINUE
10411 IF(IBOLD .GE. IBNEW) GOTO 10412
      XOUT=X(IBOLD)
      YOUT=Y(IBOLD)
      WIN=W(IBOLD)
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)
      SUMW=SUMW-WIN
      IBOLD=IBOLD+1
      GOTO 10411
10412 CONTINUE
      IF(.NOT.(CROSS.eq.1))GOTO 10431
      XOUT=X(I)
      YOUT=Y(I)
      WIN=W(I)
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)
      SUMW=SUMW-WIN
10431 CONTINUE
      IF(VAR .LE. 0.)GOTO 10451
      SMO(I)=YBAR+COV*(X(I)-XBAR)/VAR
      DOF=DOF+W(I)/SUMW+  (W(I)*(X(I)-XBAR)**2)/VAR
      GOTO 10461
10451 CONTINUE
      SMO(I)=YBAR
      DOF=DOF+W(I)/SUMW
10461 CONTINUE
      CONTINUE
      IF(.NOT.(CROSS.eq.1))GOTO 10481
      XIN=X(I)
      YIN=Y(I)
      WIN=W(I)
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW
      SUMW=SUMW+WIN
10481 CONTINUE
10371 CONTINUE
      CONTINUE
      DO 10491 I=1,N
      Y(I)=SCRAT(I)
10491 CONTINUE
      CONTINUE
      IF(CROSS .NE. 0)GOTO 10511
      I=0
10521 IF(I.GE.N-1) GOTO 10522
      I=I+1
      M0=I
10531 IF(X(I+1).GT.X(I)) GOTO 10532
      I=I+1
      IF(I .LT. N)GOTO 10531
10532 CONTINUE
      IF(I.EQ.M0)GOTO 10521
      NTIE=I-M0+1
      R=0.
      WT=0.
      DO 10541 JJ=M0,I
      J=JJ
      R=R+SMO(J)*W(J)
      WT=WT+W(J)
10541 CONTINUE
      CONTINUE
      R=R/WT
      DO 10551 J=M0,I
      SMO(J)=R
10551 CONTINUE
      CONTINUE
      GOTO 10521
10522 CONTINUE
10511 CONTINUE
      YBAR=0.0
      SUMW=0.0
      DO 10561 I=1,N
      YBAR=YBAR+W(I)*Y(I)
      SUMW=SUMW+W(I)
10561 CONTINUE
      CONTINUE
      YBAR=YBAR/SUMW
      DO 10571 I=1,N
      SMO(I)=SMO(I)-YBAR
10571 CONTINUE
      CONTINUE
      S0=YBAR
10251 CONTINUE
      CONTINUE
      RSS=0.0
      DO 10581 I=1,N
      RSS=RSS+(W(I)/SUMW)*(Y(I)-S0-SMO(I))**2
10581 CONTINUE
      CONTINUE
      
END SUBROUTINE SuperSmoother