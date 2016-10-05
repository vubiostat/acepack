C     MORTRAN 2.79 (RESERVED KEYWORD MACROS OF 09/28/81)
      SUBROUTINE RLSMO(X,Y,W,SPAN,DOF,N,SMO,RSS,SCRAT)                       22
      double precision  X(N),Y(N),W(N),SMO(N),SCRAT(N)                       23
      DOUBLE PRECISION CVRSS(6),CVSPAN(6),CVMIN, SPAN, RSS                   24
      INTEGER IDMIN                                                          25
      integer cross                                                          26
      DATA CVSPAN/0.3,0.4,0.5,0.6,0.7,1.0/                                   27
      cross=0                                                                28
      IF(span.eq.0) cross=1                                                  29
      PENAL=0.01                                                             30
      CVMIN=1E15                                                             31
      IDMIN=1                                                                33
      IF(CROSS .NE. 1)GOTO 10021                                             37
      K=1                                                                    37
      GOTO 10033                                                             37
10031 K=K+1                                                                  37
10033 IF((K).GT.(6))GOTO 10032                                               37
      CALL SMTH(X,Y,W,CVSPAN(K), DOF,N,1,SMO,S0,CVRSS(K),SCRAT)              39
      IF(CVRSS(K) .GT. CVMIN)GOTO 10051                                      40
      CVMIN=CVRSS(K)                                                         41
      IDMIN=K                                                                42
10051 CONTINUE                                                               43
      GOTO 10031                                                             44
10032 CONTINUE                                                               44
      SPAN=CVSPAN(IDMIN)                                                     45
      IF(PENAL .LE. 0.)GOTO 10071                                            46
      CVMIN=(1.+PENAL)*CVMIN                                                 48
      K=6                                                                    48
      GOTO 10083                                                             48
10081 K=K+(-1)                                                               48
10083 IF((-1)*((K)-(1)).GT.0)GOTO 10082                                      48
      IF(CVRSS(K) .GT. CVMIN)GOTO 10101                                      48
      GOTO 10082                                                             48
10101 CONTINUE                                                               49
      GOTO 10081                                                             50
10082 CONTINUE                                                               50
      SPAN=CVSPAN(K)                                                         51
10071 CONTINUE                                                               52
10021 CONTINUE                                                               53
      CALL SMTH(X,Y,W,SPAN,DOF,N,0,SMO,S0,RSS,SCRAT)                         54
      DO 10111 i=1,n                                                         54
      smo(i)=smo(i)+s0                                                       54
10111 CONTINUE                                                               55
10112 CONTINUE                                                               55
      RETURN                                                                 56
      END                                                                    57
      SUBROUTINE SMTH(X,Y,W,SPAN,DOF,N,CROSS,SMO,S0,RSS,SCRAT)               58
      double precision  X(N),Y(N),W(N),SMO(N),SCRAT(N),RSS,SPAN              59
      DOUBLE PRECISION SUMW,XBAR,YBAR,COV,VAR                                60
      INTEGER FIXEDS,CROSS,LINE                                              70
      LINE=1                                                                 83
      FIXEDS=1                                                               84
      IF(SPAN .GE. 1.0)GOTO 10131                                            84
      LINE=0                                                                 84
10131 CONTINUE                                                               86
      XBAR=X(1)                                                              86
      YBAR=Y(1)                                                              86
      COV=0.                                                                 86
      VAR=0.                                                                 86
      SUMW=W(1)                                                              88
      IF(LINE .NE. 1)GOTO 10151                                              91
      DO 10161 I=2,N                                                         91
      XIN=X(I)                                                               91
      YIN=Y(I)                                                               91
      WIN=W(I)                                                               91
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)                                    91
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)                                    91
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW                      91
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW                              91
      SUMW=SUMW+WIN                                                          91
10161 CONTINUE                                                               92
10162 CONTINUE                                                               92
      I=1                                                                    94
      GOTO 10173                                                             94
10171 I=I+1                                                                  94
10173 IF((I).GT.(N))GOTO 10172                                               94
      IF(.NOT.(CROSS.eq.1))GOTO 10191                                        94
      XOUT=X(I)                                                              94
      YOUT=Y(I)                                                              94
      WIN=W(I)                                                               94
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)                    94
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)                             94
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)                                   94
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)                                   94
      SUMW=SUMW-WIN                                                          94
10191 CONTINUE                                                               95
      IF(VAR .LE. 0.)GOTO 10211                                              97
      SMO(I)=COV*(X(I)-XBAR)/VAR                                             98
      GOTO 10221                                                             99
10211 CONTINUE                                                               99
      SMO(I)=0                                                               99
10221 CONTINUE                                                              100
10201 CONTINUE                                                              100
      IF(.NOT.(CROSS.eq.1))GOTO 10241                                       100
      XIN=X(I)                                                              100
      YIN=Y(I)                                                              100
      WIN=W(I)                                                              100
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)                                   100
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)                                   100
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW                     100
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW                             100
      SUMW=SUMW+WIN                                                         100
10241 CONTINUE                                                              101
      GOTO 10171                                                            102
10172 CONTINUE                                                              102
      S0=YBAR                                                               103
      SCRAT(1)=COV/VAR                                                      104
      DOF=1.0                                                               105
      GOTO 10251                                                            107
10151 CONTINUE                                                              111
      ITOLD=1                                                               111
      IBOLD=1                                                               111
      DOF=-1.0                                                              115
      DO 10261 I=1,N                                                        115
      SCRAT(I)=Y(I)                                                         115
10261 CONTINUE                                                              116
10262 CONTINUE                                                              116
      IF(.NOT.(cross.eq.0))GOTO 10281                                       117
      I=0                                                                   118
10291 IF(I.GE.N-1) GOTO 10292                                               118
      I=I+1                                                                 119
      M0=I                                                                  120
10301 IF(X(I+1).GT.X(I)) GOTO 10302                                         120
      I=I+1                                                                 120
      IF(I .LT. N)GOTO 10301                                                120
10302 CONTINUE                                                              121
      IF(I.EQ.M0)GOTO 10291                                                 122
      NTIE=I-M0+1                                                           123
      R=0.                                                                  123
      WT=0.                                                                 123
      DO 10311 JJ=M0,I                                                      123
      J=JJ                                                                  124
      R=R+Y(J)*W(J)                                                         124
      WT=WT+W(J)                                                            124
10311 CONTINUE                                                              124
10312 CONTINUE                                                              124
      R=R/WT                                                                125
      DO 10321 J=M0,I                                                       125
      Y(J)=R                                                                125
10321 CONTINUE                                                              126
10322 CONTINUE                                                              126
      GOTO 10291                                                            127
10292 CONTINUE                                                              127
10281 CONTINUE                                                              128
      ISPAN=N*SPAN                                                          129
      IF(.NOT.(FIXEDS.eq.1))GOTO 10341                                      129
      IS2=ISPAN/2                                                           129
      IF(IS2 .GE. 1)GOTO 10361                                              129
      IS2=1                                                                 129
10361 CONTINUE                                                              129
10341 CONTINUE                                                              135
      DO 10371 I=1,N                                                        136
      ITNEW=MIN(I+IS2,N)                                                    136
      IBNEW=MAX(I-IS2,1)                                                    137
10381 IF(ITOLD .GE. ITNEW) GOTO 10382                                       137
      ITOLD=ITOLD+1                                                         137
      XIN=X(ITOLD)                                                          137
      YIN=Y(ITOLD)                                                          137
      WIN=W(ITOLD)                                                          137
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)                                   137
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)                                   137
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW                     137
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW                             137
      SUMW=SUMW+WIN                                                         137
      GOTO 10381                                                            138
10382 CONTINUE                                                              138
10391 IF(IBOLD .LE. IBNEW) GOTO 10392                                       138
      IBOLD=IBOLD-1                                                         138
      XIN=X(IBOLD)                                                          138
      YIN=Y(IBOLD)                                                          138
      WIN=W(IBOLD)                                                          138
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)                                   138
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)                                   138
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW                     138
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW                             138
      SUMW=SUMW+WIN                                                         138
      GOTO 10391                                                            139
10392 CONTINUE                                                              139
10401 IF(ITOLD .LE. ITNEW) GOTO 10402                                       139
      XOUT=X(ITOLD)                                                         139
      YOUT=Y(ITOLD)                                                         139
      WIN=W(ITOLD)                                                          139
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)                   139
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)                            139
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)                                  139
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)                                  139
      SUMW=SUMW-WIN                                                         139
      ITOLD=ITOLD-1                                                         139
      GOTO 10401                                                            140
10402 CONTINUE                                                              140
10411 IF(IBOLD .GE. IBNEW) GOTO 10412                                       140
      XOUT=X(IBOLD)                                                         140
      YOUT=Y(IBOLD)                                                         140
      WIN=W(IBOLD)                                                          140
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)                   140
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)                            140
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)                                  140
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)                                  140
      SUMW=SUMW-WIN                                                         140
      IBOLD=IBOLD+1                                                         140
      GOTO 10411                                                            142
10412 CONTINUE                                                              142
      IF(.NOT.(CROSS.eq.1))GOTO 10431                                       142
      XOUT=X(I)                                                             142
      YOUT=Y(I)                                                             142
      WIN=W(I)                                                              142
      COV=COV-WIN*(XOUT-XBAR)*(YOUT-YBAR)*SUMW/(SUMW-WIN)                   142
      VAR=VAR-WIN*(XOUT-XBAR)**2*SUMW/(SUMW-WIN)                            142
      XBAR=(SUMW*XBAR-WIN*XOUT)/(SUMW-WIN)                                  142
      YBAR=(SUMW*YBAR-WIN*YOUT)/(SUMW-WIN)                                  142
      SUMW=SUMW-WIN                                                         142
10431 CONTINUE                                                              143
      IF(VAR .LE. 0.)GOTO 10451                                             145
      SMO(I)=YBAR+COV*(X(I)-XBAR)/VAR                                       146
      DOF=DOF+W(I)/SUMW+  (W(I)*(X(I)-XBAR)**2)/VAR                         148
      GOTO 10461                                                            149
10451 CONTINUE                                                              149
      SMO(I)=YBAR                                                           149
      DOF=DOF+W(I)/SUMW                                                     149
10461 CONTINUE                                                              150
10441 CONTINUE                                                              150
      IF(.NOT.(CROSS.eq.1))GOTO 10481                                       150
      XIN=X(I)                                                              150
      YIN=Y(I)                                                              150
      WIN=W(I)                                                              150
      XBAR=(SUMW*XBAR+XIN*WIN)/(SUMW+WIN)                                   150
      YBAR=(SUMW*YBAR+YIN*WIN)/(SUMW+WIN)                                   150
      COV=COV+WIN*(XIN-XBAR)*(YIN-YBAR)*(SUMW+WIN)/SUMW                     150
      VAR=VAR+WIN*(XIN-XBAR)**2*(SUMW+WIN)/SUMW                             150
      SUMW=SUMW+WIN                                                         150
10481 CONTINUE                                                              151
10371 CONTINUE                                                              155
10372 CONTINUE                                                              155
      DO 10491 I=1,N                                                        155
      Y(I)=SCRAT(I)                                                         155
10491 CONTINUE                                                              156
10492 CONTINUE                                                              156
      IF(CROSS .NE. 0)GOTO 10511                                            157
      I=0                                                                   158
10521 IF(I.GE.N-1) GOTO 10522                                               158
      I=I+1                                                                 159
      M0=I                                                                  160
10531 IF(X(I+1).GT.X(I)) GOTO 10532                                         160
      I=I+1                                                                 160
      IF(I .LT. N)GOTO 10531                                                160
10532 CONTINUE                                                              161
      IF(I.EQ.M0)GOTO 10521                                                 162
      NTIE=I-M0+1                                                           163
      R=0.                                                                  163
      WT=0.                                                                 163
      DO 10541 JJ=M0,I                                                      163
      J=JJ                                                                  164
      R=R+SMO(J)*W(J)                                                       164
      WT=WT+W(J)                                                            164
10541 CONTINUE                                                              164
10542 CONTINUE                                                              164
      R=R/WT                                                                165
      DO 10551 J=M0,I                                                       165
      SMO(J)=R                                                              165
10551 CONTINUE                                                              166
10552 CONTINUE                                                              166
      GOTO 10521                                                            167
10522 CONTINUE                                                              167
10511 CONTINUE                                                              168
      YBAR=0.0                                                              168
      SUMW=0.0                                                              169
      DO 10561 I=1,N                                                        169
      YBAR=YBAR+W(I)*Y(I)                                                   169
      SUMW=SUMW+W(I)                                                        169
10561 CONTINUE                                                              170
10562 CONTINUE                                                              170
      YBAR=YBAR/SUMW                                                        171
      DO 10571 I=1,N                                                        171
      SMO(I)=SMO(I)-YBAR                                                    171
10571 CONTINUE                                                              172
10572 CONTINUE                                                              172
      S0=YBAR                                                               173
10251 CONTINUE                                                              174
10141 CONTINUE                                                              174
      RSS=0.0                                                               175
      DO 10581 I=1,N                                                        175
      RSS=RSS+(W(I)/SUMW)*(Y(I)-S0-SMO(I))**2                               175
10581 CONTINUE                                                              178
10582 CONTINUE                                                              178
      RETURN                                                                178
      END                                                                   178
