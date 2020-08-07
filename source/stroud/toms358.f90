Module mconst
      Implicit None
      Public
!
      Integer, Parameter :: FP = Kind (1d0)
      Integer, Parameter :: DP = Kind (1d0)
      Integer, Parameter :: CP=kind((0d0,0d0))
      Integer, Parameter :: TIME = Kind (1d0)
      Real (FP), Parameter :: PI2 = 2._FP * Atan2 (0._FP,-1._FP)
      Character (Len=1), Parameter :: NL = NEW_LINE ("l")
Contains
      Subroutine err (test, msg, warn)
         Logical, Intent (In) :: test
         Logical, Optional, Intent (In) :: warn
         Character (Len=*), Intent (In) :: msg
!
         Logical :: warndef
         If ( .Not. test) Return
         warndef = .False.
         If (present(warn)) warndef = warn
         Print *, trim (msg)
         If ( .Not. warndef) Stop
      End Subroutine err
End Module mconst

Module msub90
      Use mconst
      Implicit None
Contains
!
!
      Subroutine CSVD (A, NU, NV, P, S, U, V)
!
! SVD factorization: U(diag(S))V^H = A, where U(MxM),S(N),V(NxN)
!
!  Discussion:
!
!    CSVD computes Singular Value Decomposition A(MxN)=U(MxM)S(MxN)V(NxN)^H
!    where returned U,V are unitary and S is null except non-negative diagonal
!    starting at S(1,1) and ^H is Hermitte complex conjugate transpose. 
!    Note: N<=M. 
!
!    Number representation parameters:
!    The original text  ETA = 1.5E-8 and TOL = 1.E-31 were replaced
!    by calls to EPSILON() and TINY().
!
!  Modified:
!
!    02 December 2017
!
!  Author:
!
!    P A Businger, G H Golub.
!    Upgrade to f90 by A.Schwarzenberg-Czerny
!
!  Reference:
!
!    Peter Businger, Gene Golub,
!    Algorithm 358:
!    Singular Value Decomposition of a Complex Matrix,
!    Communications of the ACM,
!    Volume 12, Number 10, October 1969, pages 564-565.
!    
!  Parameters:
!
!    Input:
!
!    NU = 0|N,M - The only permitted NU is 0,N or M. For P>0 required NU>=N.
!    NV = 0|N   - The only permitted values of NV are 0 and N.
!    A(M,N+P)   - complex, input/output M by N matrix. For P>0 A(N+1:N+P) = B, 
!               where AX = B is multi RHSs least squares (LSQ) system.
!    Output:
!
!    S(N)       - real non-negative singular values sorted in descending 
!                 order (and so are corresponding rows/columns of U/V). 
!                 Number of non-negligible diagonal elements of S is order of A.
!    U(M,NU)    - complex , output, the first NU columns of U. 
!    V(M,NV)    - complex , output, the first NV columns of V.
!    A(N+1:N+P) - for P>0 contain U^HB, where LSQ solution X = V(diag(1/S))U^HB
!                 and the rest of A is overwritten. X(MxP)
!
         Implicit None
         Integer, Intent (In) :: NU, NV, P
         Complex (CP), Intent (Inout) :: A (:,:)
         Complex (CP), Intent (Out) :: U (:,:), V (:,:)
         Real (FP), Intent (Out) :: S (:)
!
         Integer :: M, N, NP, N1, K, K1, I, J, L, L1
         Complex (CP) :: Q, CSWAP (Size(A, Dim=2))
         Real (FP) :: ETA=Epsilon(ETA), TOL=Tiny(TOL)/Epsilon(TOL)
         Real (FP), Dimension(Size(S)) :: B, C, T, SX, SY
         Real (FP) :: Z, W, EPS, CS, SN, F, G, H, X, Y
!
! Dimensions check
         M  = Size(A, Dim=1)
         NP = Size(A, Dim=2)
         N  = NP-P
         N1 = N + 1
         Call err( M<N .Or. P<0, "CSVD: Error: M<N .Or. P<0")
         Call err ( 0 == count ( (/0,N,M/) - NU == 0 ), "CSVD: Error: NU/=0|N|M") 
         Call err(1/=count((/0,N/)-NV==0),"CSVD: Error: NV/=0|N") 
         Call err(P>0 .And. NU<N, "CSVD: Error: P>0 .And. NU<N") 
         Call err(Size(S)/=N .Or. Size(U,Dim=1)/=M .Or. Size(U,Dim=2)/=NU &
            .Or. Size(V, Dim=1)/=NV .Or. Size(V, Dim=2)/=NV, &
            "CSVD: Error: Dimensions of SUV differ from N,MxNU,NVxNV")
! Householder reduction
         C (1) = 0._FP
         K = 1
         Do
            K1 = K + 1
! Elimination of A(I,K), I=K+1,...M
            Z = Sum(Abs (A(K:, K)) ** 2)
            B (K) = 0._FP
            If (Z > TOL) Then
               Z = Sqrt (Z)
               B (K) = Z
               W = Abs (A(K, K))
               Q = (1._FP, 0._FP)
               If (W /= 0._FP) Q = A (K, K) / W
               A (K, K) = Q * (Z+W)
               If (K /= NP) Then
                  Do J = K1, NP
                     Q = Dot_product(A(K:, K), A (K:, J))/ (Z*(Z+W))
                     A (K:, J) = A (K:, J) - Q * A (K:, K)
                  End Do
! Phase transformation
                  Q = - CONJG (A(K, K)) / Abs (A(K, K))
                  A (K, K1:NP) = Q * A (K, K1:NP)
               End If
            End If
! Elimination of A(K,J), J=K+2,...,N
            If (K == N) Exit! 140
            Z = Sum(Abs (A(K, K1:N)) ** 2)
            C (K1) = 0._FP
            If (Z > TOL) Then
               Z = Sqrt (Z)
               C (K1) = Z
               W = Abs (A(K, K1))
               Q = (1._FP, 0._FP)
               If (W /= 0._FP) Q = A (K, K1) / W
               A (K, K1) = Q * (Z+W)
               Do I = K1, M
                  Q = Dot_product(A(K, K1:N), A (I, K1:N))/ (Z*(Z+W))
                  A (I, K1:N) = A (I, K1:N) - Q * A (K, K1:N)
               End Do
! Phase transformation
               Q = - CONJG (A(K, K1)) / Abs (A(K, K1))
               A (K1:, K1) = A (K1:, K1) * Q
            End If
            K = K1
         End Do
! Tolerance for negligible elementa
         ! 140
         S = B
         T = C
         EPS = MAXVAL(S+T) * ETA
! Initialization of U and V
         U(:,:NU)=(0._FP, 0._FP)
         Forall (J = 1 : NU) U (J, J) = (1._FP, 0._FP)
         V(:,:NV)=(0._FP, 0._FP)
         Forall (J = 1 : NV) V (J, J) = (1._FP, 0._FP)
! QR diagonalization
         Do K = N, 1, -1
! Test for split
            Do ! 220 loop
               Do L = K, 1, -1
                  If (Abs(T(L)) <= EPS) Exit
                  If (Abs(S(L-1)) <= EPS) Exit
               End Do
               If (Abs(T(L)) > EPS) Then ! GoTo 290
! Cancellation of B(L)
                  CS = 0._FP
                  SN = 1._FP
                  L1 = L - 1
                  Do I = L, K
                     F = SN * T (I)
                     T (I) = CS * T (I)
                     If (Abs(F) <= EPS) Exit
                     H = S (I)
                     W = Sqrt (F*F+H*H)
                     S (I) = W
                     CS = H / W
                     SN = - F / W
                     If (NU /= 0) Then
                        SX = REAL (U(:N, L1), FP)
                        SY = REAL (U(:N, I),  FP)
                        U (:N, L1) = Cmplx (SX*CS+SY*SN, 0._FP, Kind=CP)
                        U (:N, I)  = Cmplx (SY*CS-SX*SN, 0._FP, Kind=CP)
                     End If
                     If (NP /= N) Then
                        CSWAP (N1 : NP) = A (L1, N1 : NP)
                        A (L1, N1 : NP) = A (L1, N1 : NP) * CS + &
                                          A (I,  N1 : NP) * SN
                        A (I,  N1 : NP) = A (I,  N1 : NP) * CS - &
                                          CSWAP (N1 : NP) * SN
                     End If
                  End Do
               End If
! Test for convergence
               W = S (K)
               If (L == K) Exit! GoTo 360
! Origin shift
               X = S (L)
               Y = S (K-1)
               G = T (K-1)
               H = T (K)
               F = ((Y-W)*(Y+W)+(G-H)*(G+H)) / (2._FP*H*Y)
               G = Sqrt (F*F+1._FP)
               If (F < 0._FP) G = - G
               F = ((X-W)*(X+W)+(Y/(F+G)-H)*H) / X
! QR step
               CS = 1._FP
               SN = 1._FP
               L1 = L + 1
               Do I = L1, K
                  G = T (I)
                  Y = S (I)
                  H = SN * G
                  G = CS * G
                  W = Sqrt (H*H+F*F)
                  T (I-1) = W
                  CS = F / W
                  SN = H / W
                  F = X * CS + G * SN
                  G = G * CS - X * SN
                  H = Y * SN
                  Y = Y * CS
                  If (NV /= 0) Then
                     SX = REAL (V(:, I-1),FP)
                     SY = REAL (V(:, I),  FP)
                     V (:, I-1) = Cmplx (SX*CS+SY*SN, 0._FP, Kind=CP)
                     V (:, I)   = Cmplx (SY*CS-SX*SN, 0._FP, Kind=CP)
                  End If
                  W = Sqrt (H*H+F*F)
                  S (I-1) = W
                  CS = F / W
                  SN = H / W
                  F = CS * G + SN * Y
                  X = CS * Y - SN * G
                  If (NU /= 0) Then
                     SX = REAL (U(:N, I-1),FP)
                     SY = REAL (U(:N, I),  FP)
                     U (:N, I-1) = Cmplx (SX*CS+SY*SN, 0._FP, Kind=CP)
                     U (:N, I)   = Cmplx (SY*CS-SX*SN, 0._FP, Kind=CP)
                  End If
                  If (N /= NP) Then
                     CSWAP  (N1 : NP) = A (I-1, N1 : NP)
                     A (I-1, N1 : NP) = A (I-1, N1 : NP) * CS + &
                                        A (I,   N1 : NP) * SN
                     A (I,   N1 : NP) = A (I,   N1 : NP) * CS - &
                                        CSWAP  (N1 : NP) * SN
                  End If
               End Do
               T (L) = 0._FP
               T (K) = F
               S (K) = X
            End Do ! goto 220
! Convergence
            If (W < 0._FP) Then ! 360
               S (K) = - W      ! Returns positive singular values
               If (NV /= 0) V(:,K) = - V(:,K)
            End If
         End Do ! 380
! Sort singular values
         Do K = 1, N
            J = Maxloc( S(K:), Dim = 1) + K - 1
            G = S(J)
            If (J /= K) Then
               S (J) = S (K)
               S (K) = G
               If (NV /= 0) Then ! Interchange V(1:N,J) and V(1:N,K).
                  CSWAP (:N) = V (:, J)
                  V (:, J)   = V (:, K)
                  V (:, K) = CSWAP (:N) 
               End If
               If (NU /= 0) Then ! Interchange U(1:N,J) and U(1:N,K).
                  CSWAP (:N) = U (:N, J)
                  U (:N, J)  = U (:N, K)
                  U (:N, K) = CSWAP (:N) 
               End If            ! Interchange A(J,N1:NP) and A(K,N1:NP).
               If (N /= NP) Then
                  CSWAP (N1 : NP) =  A (J, N1 : NP)
                  A  (J, N1 : NP) =  A (K, N1 : NP)
                  A  (K, N1 : NP) = CSWAP (N1 : NP) 
               End If
            End If
         End Do
! Back transformation
         If (NU /= 0) Then
            Do K = N, 1, -1
               If (B(K) /= 0._FP) Then
                  Q = - A (K, K) / Abs (A(K, K))
                  U (K, :NU) = Q * U (K, :NU)
                  Do J = 1, NU
                     Q = Dot_Product(A(K:, K), U(K:, J))/ &
                         (Abs(A(K, K))*B(K))
                     U(K:, J) = U(K:, J) - Q * A (K:, K)
                  End Do
               End If
            End Do
         End If
         If (NV /= 0) Then
            Do K = N - 1, 1, -1
               K1 = K + 1
               If (C(K1) /= 0._FP) Then
                  Q = - CONJG (A(K, K1)) / Abs (A(K, K1))
                  V (K1, :NV) = Q * V (K1, :NV)
                  Do J = 1, NV
                     Q = Sum(A (K, K1:N) * V (K1:, J))/ &
                         (Abs(A(K, K1))*C(K1))
                     V (K1:, J) = V (K1:, J) - Q * CONJG (A(K, K1:N))
                  End Do
               End If
            End Do
         End If
      End Subroutine CSVD
!
!
End Module msub90
!
