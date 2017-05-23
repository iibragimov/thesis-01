        !COMPILER-GENERATED INTERFACE MODULE: Wed Feb 24 23:24:08 2016
        MODULE SOLVE_QUETION__genmod
          INTERFACE 
            SUBROUTINE SOLVE_QUETION(F1,F2,DF1_DX,DF1_DY,DF2_DX,DF2_DY, &
     &X0,Y0)
              REAL(KIND=8) :: F1
              EXTERNAL F1
              REAL(KIND=8) :: F2
              EXTERNAL F2
              REAL(KIND=8) :: DF1_DX
              EXTERNAL DF1_DX
              REAL(KIND=8) :: DF1_DY
              EXTERNAL DF1_DY
              REAL(KIND=8) :: DF2_DX
              EXTERNAL DF2_DX
              REAL(KIND=8) :: DF2_DY
              EXTERNAL DF2_DY
              REAL(KIND=8) :: X0
              REAL(KIND=8) :: Y0
            END SUBROUTINE SOLVE_QUETION
          END INTERFACE 
        END MODULE SOLVE_QUETION__genmod
