MODULE MODULE_CREATE_SURFACE
    
    TYPE SURFACE
        INTEGER :: NG     ! NUMBER OF GRIDS
        INTEGER :: NG1    ! NUMBER OF GRIDS - SPAN DIRECTION
        INTEGER :: NG2    ! NUMBER OF GRIDS - CHORD DIRECTION
        INTEGER :: NP     ! NUMBER OF PANELS
        REAL*8, ALLOCATABLE :: GRIDS(:,:)      ! NUMBER OF GRIDS | X,Y,Z,SPAN,CHORD
        INTEGER, ALLOCATABLE :: CN(:,:)        ! NUMBER OF PANELS - CONNECTIVITY
        REAL*8, ALLOCATABLE :: PANELS(:,:)     ! NUMBER OF PANELS - X,Y,Z,AXY,AXZ,AYZ
    END TYPE SURFACE
    
    
    TYPE CONFIG
        REAL*8, ALLOCATABLE :: LE(:,:)  ! NUMBER OF POINTS DEFINING THE LEADING EDGE | X,Y,Z,CHORD,INCIDENCE
        REAL*8, ALLOCATABLE :: TE(:,:)
        INTEGER :: DIRS, MDIR
        INTEGER :: N      ! SPAN STRIPS
        INTEGER :: M      ! CHORD STRIPS
        INTEGER :: COPT, SOPT
        INTEGER :: POPT
        REAL*8, ALLOCATABLE :: CDEF(:), SDEF(:)
        INTEGER :: NLE
        REAL*8 :: DX, DY, DZ
        TYPE(AIRFOIL), ALLOCATABLE :: PROFILES(:)
    END TYPE CONFIG
    
    
    TYPE AIRFOIL
        REAL*8, ALLOCATABLE :: UPSUF(:,:), LOSUF(:,:)
        REAL*8, ALLOCATABLE :: UPSUF2(:,:), LOSUF2(:,:)
        INTEGER :: NUP, NLO
        CHARACTER(LEN=30) :: FNAME
    END TYPE
    
    
    CONTAINS
    
    !============================================================
    SUBROUTINE WRITE_VTK(S)
        TYPE(SURFACE), INTENT(IN) :: S
        INTEGER :: I
        
        OPEN(UNIT=40, FILE='surface.vtk', STATUS='UNKNOWN')
        
        WRITE(40, '("# vtk DataFile Version 2.0")')
        WRITE(40, '("surface.vtk")')
        WRITE(40, '("ASCII")')
        WRITE(40, '("DATASET UNSTRUCTURED_GRID")')
        !WRITE(40, *)
        
        WRITE(40, '("POINTS",I8," float")')S%NG
        
        DO I=1, S%NG
            WRITE(40,'(3F12.6)')S%GRIDS(I,1), S%GRIDS(I,2), S%GRIDS(I,3)
        END DO
        
        WRITE(40,'("CELLS ",I8,I8)')S%NP, S%NP*5
        
        DO I=1, S%NP
            WRITE(40,'(I2,4I8)')4, S%CN(I,1)-1, S%CN(I,2)-1, S%CN(I,3)-1, S%CN(I,4)-1
        END DO

        !WRITE(40, *)
        
        WRITE(40,'("CELL_TYPES ",I8)')S%NP
        
        DO I=1, S%NP
            WRITE(40,'(I1)')9
        END DO
        
        
        CLOSE(40)
    
    END SUBROUTINE
    
    !============================================================    
    SUBROUTINE DEF_SURFACE(C, S)
    
        TYPE(CONFIG), INTENT(INOUT) :: C
        TYPE(SURFACE), INTENT(INOUT) :: S
        
        REAL*8 :: PT(3,2), LCHORD
        INTEGER :: I, J, K, L
        
        S%NP = C%N * C%M
        S%NG1 = C%N + 1
        S%NG2 = C%M + 1
        S%NG = S%NG1 * S%NG2 
        
        ALLOCATE(S%PANELS(S%NP,6))
        ALLOCATE(S%CN(S%NP,4))
        ALLOCATE(S%GRIDS(S%NG ,5))
        
        OPEN(UNIT=11, FILE='grids_log.txt', STATUS='UNKNOWN')
        OPEN(UNIT=111, FILE='grids.txt', STATUS='UNKNOWN')
        OPEN(UNIT=12, FILE='connectivity.txt', STATUS='UNKNOWN') 
        OPEN(UNIT=13, FILE='panels.txt', STATUS='UNKNOWN')
        OPEN(UNIT=21, FILE='log.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        
        WRITE(21,*)
        WRITE(21,'("Grids spanwise:  ",I8)')S%NG1
        WRITE(21,'("Grids clockwise: ",I8)')S%NG2
        WRITE(21,'("Panels:          ",I8)')S%NP
        WRITE(21,'("  Grid    Span    Chord       xle         yle         zle         xte         yte         zte         dx")')
        WRITE(11,'("  Grid    Span    Chord       x[m]        y[m]        z[m]  ")')
        WRITE(111,'(I8)')S%NG
        WRITE(12, '(I8)')S%NP
        
        K = 0
        
        ! Grids
        DO I=1, S%NG2   ! chordwise
            DO J=1, S%NG1  ! spanwise
                K = K + 1
                S%GRIDS(K, 4) = J
                S%GRIDS(K, 5) = I
                IF (C%DIRS == 2) THEN
                    IF (C%SOPT == 0) THEN
                        S%GRIDS(K, 2) = C%LE(1,2) + (J - 1) * C%DY
                    ELSE IF (C%SOPT == 1) THEN
                        S%GRIDS(K, 2) = C%LE(1,2) + (C%LE(C%NLE, 2) - C%LE(1, 2)) * C%SDEF(J)
                    END IF
                    PT(2,1) = S%GRIDS(K, 2)
                    PT(2,2) = S%GRIDS(K, 2)
                    
                    ! Point at the Leading Edge
                    CALL LININTERP(C%LE(:,2),C%LE(:,1),C%NLE,PT(2,1),PT(1,1))
                    CALL LININTERP(C%LE(:,2),C%LE(:,3),C%NLE,PT(2,1),PT(3,1))
                    ! Point at the Trailing Edge
                    CALL LININTERP(C%TE(:,2),C%TE(:,1),C%NLE,PT(2,2),PT(1,2))
                    CALL LININTERP(C%TE(:,2),C%TE(:,3),C%NLE,PT(2,2),PT(3,2))
                    IF (C%COPT == 0) THEN
                        C%DX = (I - 1) * (PT(1,2) - PT(1,1)) / S%NG2
                        S%GRIDS(K, 1) = PT(1,1) + C%DX
                    ELSE IF (C%COPT == 1) THEN
                        S%GRIDS(K, 1) = PT(1,1) + (PT(1,2) - PT(1,1)) * C%CDEF(I)
                    END IF
                    CALL LININTERP(PT(1,:),PT(3,:),2,S%GRIDS(K, 1),S%GRIDS(K, 3))
                ELSE IF (C%DIRS == 3) THEN
                    IF (C%SOPT == 0) THEN
                        S%GRIDS(K, 3) = C%LE(1,3) + (J - 1) * C%DZ
                    ELSE IF (C%SOPT == 1) THEN
                        S%GRIDS(K, 3) = C%LE(1,3) + (C%LE(NLE, 3) - C%LE(1, 3)) * C%SDEF(J)
                    END IF                        
                    PT(3,1) = S%GRIDS(K, 3)
                    PT(3,2) = S%GRIDS(K, 3)
                    ! Point at the Leading Edge
                    CALL LININTERP(C%LE(:,3),C%LE(:,1),C%NLE,PT(3,1),PT(1,1))
                    CALL LININTERP(C%LE(:,3),C%LE(:,2),C%NLE,PT(3,1),PT(2,1))
                    ! Point at the Trailing Edge
                    CALL LININTERP(C%TE(:,3),C%TE(:,1),C%NLE,PT(3,2),PT(1,2))
                    CALL LININTERP(C%TE(:,3),C%TE(:,2),C%NLE,PT(3,2),PT(2,2))
                    C%DX = (I - 1) * (PT(1,2) - PT(1,1)) / S%NG2
                    S%GRIDS(K, 1) = PT(1,1) + C%DX
                    CALL LININTERP(PT(1,:),PT(2,:),2,S%GRIDS(K, 1),S%GRIDS(K, 2))
                END IF
                WRITE(11,'(3I8,3F12.4)')K, J, I, S%GRIDS(K, 1), S%GRIDS(K, 2), S%GRIDS(K, 3)
                WRITE(111,'(I8,3F12.4)')K, S%GRIDS(K, 1), S%GRIDS(K, 2), S%GRIDS(K, 3)
                WRITE(21,'(3I8,7F12.4)')K, J, I, PT(1,1), PT(2,1), PT(3,1), PT(1,2), PT(2,2), PT(3,2), C%DX
            END DO
        END DO
        
        K = 0 ! Spanwise
        J = 1 ! Chordwise
        !Connectivity
        DO I = 1, S%NP
            IF (K == C%N) THEN !Span strip
                K = 1
                J = J + 1
            ElSE
                K = K + 1
            END IF
            IF (C%MDIR == 1) THEN
                S%CN(I,1) = I + (J - 1)
                S%CN(I,2) = S%NG1 + I + (J - 1)
                S%CN(I,3) = S%NG1 + I + 1 + (J - 1)
                S%CN(I,4) = I + 1 + (J - 1)
            ELSE IF (C%MDIR == 0) THEN
                S%CN(I,1) = I + (J - 1)
                S%CN(I,4) = S%NG1 + I + (J - 1)
                S%CN(I,3) = S%NG1 + I + 1 + (J - 1)
                S%CN(I,2) = I + 1 + (J - 1)
            END IF
            WRITE(12, '(4I8)')S%CN(I,1),S%CN(I,2),S%CN(I,3),S%CN(I,4)
        END DO
        
        
        !Panels
        DO I = 1, S%NP
            S%PANELS(I,1) = 0.25*(S%GRIDS(S%CN(I,1),1) + S%GRIDS(S%CN(I,2),1) + S%GRIDS(S%CN(I,3),1) + S%GRIDS(S%CN(I,4),1))
            S%PANELS(I,2) = 0.25*(S%GRIDS(S%CN(I,1),2) + S%GRIDS(S%CN(I,2),2) + S%GRIDS(S%CN(I,3),2) + S%GRIDS(S%CN(I,4),2))
            S%PANELS(I,3) = 0.25*(S%GRIDS(S%CN(I,1),3) + S%GRIDS(S%CN(I,2),3) + S%GRIDS(S%CN(I,3),3) + S%GRIDS(S%CN(I,4),3))
            CALL AREA_QUAD(S%GRIDS(S%CN(I,1),1), S%GRIDS(S%CN(I,2),1), S%GRIDS(S%CN(I,3),1), S%GRIDS(S%CN(I,4),1),  &
                           S%GRIDS(S%CN(I,1),2), S%GRIDS(S%CN(I,2),2), S%GRIDS(S%CN(I,3),2), S%GRIDS(S%CN(I,4),2),  &
                           S%GRIDS(S%CN(I,1),3), S%GRIDS(S%CN(I,2),3), S%GRIDS(S%CN(I,3),3), S%GRIDS(S%CN(I,4),3),  &
                           S%PANELS(I,4), S%PANELS(I,5), S%PANELS(I,6))
            
            WRITE(13,'(I8,3F12.4,3F12.8)')I, S%PANELS(I,1:6)
            
        END DO
        
        
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
    
    END SUBROUTINE
    !============================================================
    
    

    !============================================================
    SUBROUTINE SET_CONFIG(C)
    
        TYPE(CONFIG), INTENT(INOUT) :: C
        INTEGER :: I, J
        REAL*8 :: DEG2RAD
        
        DEG2RAD = ACOS(-1.0) / 180
        
        OPEN(UNIT=10, FILE='config.txt', STATUS='OLD')
        OPEN(UNIT=21, FILE='log.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        
        READ(10,*)C%MDIR
        READ(10,*)C%DIRS
        READ(10,*)C%N
        READ(10,*)C%M
        READ(10,*)C%POPT
        READ(10,*)C%COPT
        IF (C%COPT == 1) THEN
            ALLOCATE(C%CDEF(C%M + 1))  ! GRIDS POSITIONS CHORDWISE
            READ(10, *)(C%CDEF(I), I=1,SIZE(C%CDEF))
        END IF
        READ(10,*)C%SOPT
        IF (C%SOPT == 1) THEN
            ALLOCATE(C%SDEF(C%N + 1))  ! GRIDS POSITIONS SPANWISE
            READ(10, *)(C%SDEF(I), I=1,SIZE(C%SDEF))
        END IF        
        READ(10,*)C%NLE
        
        IF (C%POPT == 1) then
            ALLOCATE(C%PROFILES(C%NLE))
        END IF
        ALLOCATE(C%LE(C%NLE, 5))
        ALLOCATE(C%TE(C%NLE, 3))
        
        WRITE(21,'("Leading and Trailing Edges:")')
        WRITE(21,'("      xle          yle          zle          xte          yte          zte")')
        
        DO I=1,C%NLE
            
            IF (C%POPT == 1) THEN
                READ(10,*)(C%LE(I,J),J=1,5), C%PROFILES(I)%FNAME
                CALL READ_AIRFOIL(C%PROFILES(I))
            ELSE
                READ(10,*)(C%LE(I,J),J=1,5)
            END IF
            
            C%TE(I,1) = C%LE(I,1) + C%LE(I,4)
            IF (C%DIRS == 2) THEN
                C%TE(I,2) = C%LE(I,2)
                IF (C%LE(I,5) .NE. 0) THEN
                    C%TE(I,3) = C%LE(I,3) - C%LE(I,4) * TAN(C%LE(I,5) * DEG2RAD)
                ELSE
                    C%TE(I,3) = C%LE(I,3)
                END IF
            ELSE IF (C%DIRS == 3) THEN
                C%TE(I,3) = C%LE(I,3)
                IF (C%LE(I,5) .NE. 0) THEN
                    C%TE(I,2) = C%LE(I,2) - C%LE(I,4) * TAN(C%LE(I,5) * DEG2RAD) 
                ELSE
                    C%TE(I,2) = C%LE(I,2)
                END IF
            END IF
            
            WRITE(21, '(6(F12.4,1X))')C%LE(I,1:3),C%TE(I,1:3)
            
        END DO
        
        IF (C%DIRS == 2) THEN
            C%DY = (C%LE(C%NLE,2) - C%LE(1,2)) / C%N
        ELSE IF (C%DIRS == 3) THEN
            C%DZ = (C%LE(C%NLE,3) - C%LE(1,3)) / C%N
        END IF
            
        
        
        CLOSE(10)
        CLOSE(21)
    
    END SUBROUTINE
    !============================================================

    
    
    !===================================================
    !
    ! LINEAR INTERPOLATION
    !
    !==================================================   
    SUBROUTINE LININTERP(X,Y,NP,XEX,YS)
         
         INTEGER,INTENT(IN) :: NP
         REAL*8,INTENT(IN) :: X(NP), Y(NP), XEX
         REAL*8,INTENT(OUT) :: YS


         IF(X(2) > X(1)) THEN                                                                                      
            IF(XEX < X(1)) THEN
              !YS = Y(1)
              YS = (Y(2)-Y(1))*(XEX-X(1))/(X(2)-X(1)) + Y(1)                                                         
            ELSE                                                                                                     
               IF(XEX >= X(NP)) THEN 
                  !YS = Y(NP)
                  YS = (Y(NP)-Y(NP-1))*(XEX-X(NP))/(X(NP)-X(NP-1))+Y(NP)                                               
               ELSE                                                                                                   
                  DO I = 2, NP                                                                                         
                     IF(XEX <= X(I)) THEN                                                                               
                        YS = (Y(I)-Y(I-1))*(XEX-X(I-1))/(X(I)-X(I-1))+Y(I-1)                                             
                        EXIT                                                                    
                     ENDIF                                                                                              
                  ENDDO                                                                                                
               ENDIF                                                                                                  
            ENDIF                                                                                                    
         ELSE                                                                                                       
            IF(XEX > X(1)) THEN
                !YS = Y(1)
                YS  = (Y(2)-Y(1))*(XEX-X(1))/(X(2)-X(1)) + Y(1)                                                        
            ELSE                                                                                                     
                IF(XEX <= X(NP)) THEN
                   !YS = Y(NP)
                   YS  = (Y(NP)-Y(NP-1))*(XEX-X(NP))/(X(NP)-X(NP-1))+Y(NP)                                              
                ELSE                                                                                                   
                   DO I = 2, NP                                                                                         
                       IF(XEX >= X(I)) THEN                                                                               
                           YS = (Y(I)-Y(I-1))*(XEX-X(I-1))/(X(I)-X(I-1))+Y(I-1)                                             
                           EXIT                                                                    
                       ENDIF                                                                                              
                   ENDDO                                                                                                
                ENDIF                                                                                                  
            ENDIF                                                                                                    
         ENDIF  

    END SUBROUTINE
    !==================================================
    
    !==================================================
    SUBROUTINE AREA_QUAD(X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, AXY, AXZ, AYZ)
    
        REAL*8, INTENT(IN) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4
        REAL*8, INTENT(INOUT) :: AXY, AXZ, AYZ
        REAL*8 :: A1, A2, A3, B1, B2, B3, VX1, VY1, VZ1
        
    !     THIS ROUTINE CALCULATE THE PANEL AREAS

    !     PRIMARY DIAGONAL (1-3)
          A1 = X3 - X1
          A2 = Y3 - Y1
          A3 = Z3 - Z1

    !      SECUNDARY DIAGONAL (2-4)
          B1 = X4 - X2
          B2 = Y4 - Y2
          B3 = Z4 - Z2

    !     VECTORIAL PRODUCT (A X B)
          VX1 = (A2 * B3 - A3 * B2)
          VY1 = (A3 * B1 - A1 * B3)
          VZ1 = (A1 * B2 - A2 * B1)
          
    !     PROJECTED AREAS AT EACH PLANE
          AXY = 0.5 * VZ1
          AXZ = 0.5 * VY1
          AYZ = 0.5 * VX1
          
    END SUBROUTINE
    !==================================================
    
    
    
    !==================================================
    SUBROUTINE READ_AIRFOIL(AF)
    
        TYPE(AIRFOIL), INTENT(INOUT) :: AF
        INTEGER :: I
        
        OPEN(UNIT=77, FILE=FNAME, STATUS='OLD')
        
        READ(77, *) AF%NUP
        ALLOCATE(AF%UPSUF(AF%NUP, 2))
        ALLOCATE(AF%UPSUF2(AF%NUP, 3))
            
        DO I = 1, AF%NUP
           READ(77, *)AF%UPSUF(I, 1), AF%UPSUF(I, 2)
        END DO

        READ(77, *) AF%NUP
        ALLOCATE(AF%UPSUF(AF%NLO, 2))
        ALLOCATE(AF%UPSUF2(AF%NLO, 3))
            
        DO I = 1, AF%NUP
           READ(77, *)AF%LOSUF(I, 1), AF%LOSUF(I, 2)
        END DO
        
        CLOSE(77)
    
    END SUBROUTINE
    
    
END MODULE