      PROGRAM TRUSS2D
C=============================================================
C  Trelica 2D - Fortran 77 (didatico)
C  GDL por no: (ux, uy)
C  Elemento: barra 2 nos (tracao/compressao)
C  CC: deslocamentos prescritos via:
C      METHOD=1: penalidade
C      METHOD=2: eliminacao de linhas/colunas
C=============================================================
      IMPLICIT NONE
      INTEGER GDL
      INTEGER MAXN, MAXE, MAXM, MAXBC, MAXL
      PARAMETER (MAXN=200, MAXE=400, MAXM=50, MAXBC=400, MAXL=400)
      
      INTEGER MAX2N
      PARAMETER (MAX2N = 2*MAXN) 

      INTEGER NN, NE, NM, NBC, NLOAD, METHOD
      INTEGER I, J, K, E, ID, JD, M
      INTEGER N1, N2, MATID
      INTEGER NODE
      INTEGER BCNODE(MAXBC), BCDOF(MAXBC)
      DOUBLE PRECISION BCVAL(MAXBC)
      INTEGER LDNODE(MAXL), LDDOF(MAXL)
      DOUBLE PRECISION LDVAL(MAXL)

      DOUBLE PRECISION X(MAXN), Y(MAXN)
      INTEGER E1(MAXE), E2(MAXE), EMAT(MAXE)
      DOUBLE PRECISION EA(MAXM), AA(MAXM)

      INTEGER NDOF, NFREE
      INTEGER G1, G2, G3, G4
      DOUBLE PRECISION KGL(4,4), KE(4,4)
      DOUBLE PRECISION C, S, L, EMod, Area, kfac

      DOUBLE PRECISION KMAT(2*MAXN,2*MAXN)
      DOUBLE PRECISION F(2*MAXN), U(2*MAXN), R(2*MAXN)
      DOUBLE PRECISION PENA
      
      DOUBLE PRECISION UE(4), DU, NAX, SIG

      DOUBLE PRECISION KMAT_ORIG(2*MAXN,2*MAXN), F_ORIG(2*MAXN)
      INTEGER IS_PRESCRIBED(2*MAXN)
      INTEGER FREE_GDLS(2*MAXN)
      DOUBLE PRECISION KRED(2*MAXN,2*MAXN), FRED(2*MAXN), URED(2*MAXN)

      CHARACTER*120 INFILE, OUTFILE
      INTEGER OUTUNIT
C-------------------------------------------------------------

      INFILE = 'truss2d02.dat'
      OUTFILE = 'truss2d02-eliminacao.out'
      OUTUNIT = 20

      OPEN(10,FILE=INFILE,STATUS='OLD')
      OPEN(OUTUNIT,FILE=OUTFILE,STATUS='UNKNOWN')

C----- Leitura basica
      READ(10,*) NN, NE, NM, NBC, NLOAD
      READ(10,*) METHOD
      NDOF = 2*NN

C----- Nos
      DO I=1,NN
         READ(10,*) NODE, X(NODE), Y(NODE)
      END DO

C----- Materiais (id, E, A)
      DO I=1,NM
         READ(10,*) M, EA(M), AA(M)
      END DO

C----- Elementos (id, no1, no2, matid)
      DO E=1,NE
         READ(10,*) ID, E1(ID), E2(ID), EMAT(ID)
      END DO

C----- Condicoes de contorno (no, dof(1=ux,2=uy), valor)
      DO I=1,NBC
         READ(10,*) BCNODE(I), BCDOF(I), BCVAL(I)
      END DO

C----- Cargas nodais (no, dof(1=ux,2=uy), valor)
      DO I=1,NLOAD
         READ(10,*) LDNODE(I), LDDOF(I), LDVAL(I)
      END DO

      CLOSE(10)

C----- Zera matrizes e vetores
      DO I=1,NDOF
         F(I)=0.0D0
         U(I)=0.0D0
         R(I)=0.0D0
         IS_PRESCRIBED(I) = 0
         DO J=1,NDOF
            KMAT(I,J)=0.0D0
         END DO
      END DO

C----- Marca gdls prescritos
      DO I=1,NBC
         G1 = GDL(BCNODE(I), BCDOF(I))
         IS_PRESCRIBED(G1) = 1
      END DO

C----- Monta vetor de cargas F
      DO I=1,NLOAD
         G1 = GDL(LDNODE(I), LDDOF(I))
         F(G1) = F(G1) + LDVAL(I)
      END DO

C=============================================================
C  Montagem global
C=============================================================
      DO E=1,NE
         N1 = E1(E)
         N2 = E2(E)
         MATID = EMAT(E)
         EMod = EA(MATID)
         Area = AA(MATID)

         L = DSQRT( (X(N2)-X(N1))**2 + (Y(N2)-Y(N1))**2 )
         C = (X(N2)-X(N1))/L
         S = (Y(N2)-Y(N1))/L

         kfac = EMod*Area/L

C----- Matriz local em coordenadas globais (4x4)
         KE(1,1) =  kfac*C*C
         KE(1,2) =  kfac*C*S
         KE(1,3) = -kfac*C*C
         KE(1,4) = -kfac*C*S

         KE(2,1) =  kfac*C*S
         KE(2,2) =  kfac*S*S
         KE(2,3) = -kfac*C*S
         KE(2,4) = -kfac*S*S

         KE(3,1) = -kfac*C*C
         KE(3,2) = -kfac*C*S
         KE(3,3) =  kfac*C*C
         KE(3,4) =  kfac*C*S

         KE(4,1) = -kfac*C*S
         KE(4,2) = -kfac*S*S
         KE(4,3) =  kfac*C*S
         KE(4,4) =  kfac*S*S

C----- Mapeamento gdl
         G1 = GDL(N1,1)
         G2 = GDL(N1,2)
         G3 = GDL(N2,1)
         G4 = GDL(N2,2)

C----- Montagem
         CALL ASSEMBLE4(KMAT, MAX2N, KE, G1,G2,G3,G4)
      END DO

C----- Armazena K e F originais
      DO I=1,NDOF
         F_ORIG(I) = F(I)
         DO J=1,NDOF
            KMAT_ORIG(I,J) = KMAT(I,J)
         END DO
      END DO

C=============================================================
C  Resolve de acordo com METHOD
C=============================================================
      IF (METHOD.EQ.1) THEN
C----- Metodo 1: Penalidade
         CALL RESOLVE_PENALTY(NDOF, MAX2N, KMAT, F, U, 
     .                        KMAT_ORIG, F_ORIG, NBC, 
     .                        BCNODE, BCDOF, BCVAL, R)
      ELSE IF (METHOD.EQ.2) THEN
C----- Metodo 2: Eliminacao de linhas/colunas
         CALL RESOLVE_ELIMINATION(NDOF, MAX2N, KMAT, F, U, 
     .                            KMAT_ORIG, F_ORIG, NBC, 
     .                            BCNODE, BCDOF, BCVAL, R,
     .                            IS_PRESCRIBED, FREE_GDLS, 
     .                            NFREE, KRED, FRED, URED)
      ELSE
         WRITE(*,*) 'Erro: METHOD deve ser 1 ou 2. Programa abortado.'
         STOP
      END IF

C=============================================================
C  Saida (escrita em arquivo OUTFILE via OUTUNIT)
C=============================================================
      WRITE(OUTUNIT,*) '========================================'
      WRITE(OUTUNIT,*) 'TRELICA 2D - RESULTADOS'
      WRITE(OUTUNIT,*) 'Arquivo de entrada: ', INFILE
      WRITE(OUTUNIT,*) 'NN, NE = ', NN, NE
      IF (METHOD.EQ.1) THEN
         WRITE(OUTUNIT,*) 'Metodo: PENALIDADE'
      ELSE
         WRITE(OUTUNIT,*) 'Metodo: ELIMINACAO DE LINHAS/COLUNAS'
      END IF
      WRITE(OUTUNIT,*) '----------------------------------------'
      WRITE(OUTUNIT,*) 'DESLOCAMENTOS (por no):'
      DO I=1,NN
         WRITE(OUTUNIT,'(A,I4,A,1PE12.4,A,1PE12.4)') 'No ',I,': ux=',
     .U(GDL(I,1)),'  uy=',U(GDL(I,2))
      END DO

      WRITE(OUTUNIT,*) '----------------------------------------'
      WRITE(OUTUNIT,*) 'REACOES (apenas gdl prescritos):'
      DO I=1,NBC
         WRITE(OUTUNIT,'(A,I4,A,I2,A,1PE12.4)') 'No ',BCNODE(I),
     .   ' dof ',BCDOF(I),'  R=',R(GDL(BCNODE(I),BCDOF(I)))
      END DO

      WRITE(OUTUNIT,*) '----------------------------------------'
      WRITE(OUTUNIT,*) 'FORCA AXIAL E TENSAO POR ELEMENTO:'
      DO E=1,NE
         N1 = E1(E)
         N2 = E2(E)
         MATID = EMAT(E)
         EMod = EA(MATID)
         Area = AA(MATID)

         L = DSQRT( (X(N2)-X(N1))**2 + (Y(N2)-Y(N1))**2 )
         C = (X(N2)-X(N1))/L
         S = (Y(N2)-Y(N1))/L

C----- deslocamentos do elemento
         Ue(1) = U(GDL(N1,1))
         Ue(2) = U(GDL(N1,2))
         Ue(3) = U(GDL(N2,1))
         Ue(4) = U(GDL(N2,2))

C----- alongamento axial (projecao)
         du  = (Ue(3)-Ue(1))*C + (Ue(4)-Ue(2))*S
         Nax = EMod*Area/L * du
         Sig = Nax/Area

         WRITE(OUTUNIT,'(A,I4,A,1PE12.4,A,1PE12.4)') 'Elem ',E,
     .   ': N=',Nax,'  sigma=',Sig
      END DO

      CLOSE(OUTUNIT)

      STOP
      END

C=============================================================
      INTEGER FUNCTION GDL(NODE, DOF)
      IMPLICIT NONE
      INTEGER NODE, DOF
      GDL = 2*(NODE-1) + DOF
      RETURN
      END

C=============================================================
      SUBROUTINE RESOLVE_PENALTY(NDOF, MAX2N, KMAT, F, U, 
     .                           KMAT_ORIG, F_ORIG, NBC, 
     .                           BCNODE, BCDOF, BCVAL, R)
      IMPLICIT NONE
      INTEGER NDOF, MAX2N, NBC
      INTEGER BCNODE(NBC), BCDOF(NBC)
      DOUBLE PRECISION BCVAL(NBC)
      DOUBLE PRECISION KMAT(MAX2N,MAX2N), F(MAX2N), U(MAX2N)
      DOUBLE PRECISION KMAT_ORIG(MAX2N,MAX2N), F_ORIG(MAX2N)
      DOUBLE PRECISION R(MAX2N)
      INTEGER GDL
      INTEGER I, J, G1
      DOUBLE PRECISION PENA

C----- Penalidade tipica: 1e12 * max(diagonal)
      PENA = 0.0D0
      DO I=1,NDOF
         IF (DABS(KMAT(I,I)).GT.PENA) PENA = DABS(KMAT(I,I))
      END DO
      IF (PENA.EQ.0.0D0) PENA = 1.0D0
      PENA = 1.0D12*PENA

      DO I=1,NBC
         G1 = GDL(BCNODE(I), BCDOF(I))
         KMAT(G1,G1) = KMAT(G1,G1) + PENA
         F(G1)       = F(G1)       + PENA*BCVAL(I)
      END DO

C----- Resolve: KMAT * U = F
      CALL GAUSS(NDOF, MAX2N, KMAT, F, U)

C----- Reacoes corretas: R = K_original * U - F_original
      DO I=1,NDOF
         R(I) = 0.0D0
         DO J=1,NDOF
            R(I) = R(I) + KMAT_ORIG(I,J)*U(J)
         END DO
         R(I) = R(I) - F_ORIG(I)
      END DO

      RETURN
      END

C=============================================================
      SUBROUTINE RESOLVE_ELIMINATION(NDOF, MAX2N, KMAT, F, U, 
     .                               KMAT_ORIG, F_ORIG, NBC, 
     .                               BCNODE, BCDOF, BCVAL, R,
     .                               IS_PRESCRIBED, FREE_GDLS, 
     .                               NFREE, KRED, FRED, URED)
      IMPLICIT NONE
      INTEGER NDOF, MAX2N, NBC, NFREE
      INTEGER BCNODE(NBC), BCDOF(NBC)
      DOUBLE PRECISION BCVAL(NBC)
      DOUBLE PRECISION KMAT(MAX2N,MAX2N), F(MAX2N), U(MAX2N)
      DOUBLE PRECISION KMAT_ORIG(MAX2N,MAX2N), F_ORIG(MAX2N)
      DOUBLE PRECISION R(MAX2N)
      INTEGER IS_PRESCRIBED(NDOF), FREE_GDLS(NDOF)
      DOUBLE PRECISION KRED(MAX2N,MAX2N), FRED(MAX2N), URED(MAX2N)
      INTEGER GDL
      INTEGER I, J, II, JJ, G1, G2
      DOUBLE PRECISION IRED, JRED

C----- Identifica gdls livres
      NFREE = 0
      DO I=1,NDOF
         IF (IS_PRESCRIBED(I).EQ.0) THEN
            NFREE = NFREE + 1
            FREE_GDLS(NFREE) = I
         END IF
      END DO

C----- Constroi sistema reduzido (apenas gdls livres)
      DO II=1,NFREE
         I = FREE_GDLS(II)
         FRED(II) = F(I)
         DO JJ=1,NFREE
            J = FREE_GDLS(JJ)
            KRED(II,JJ) = KMAT(I,J)
         END DO
      END DO

C----- Resolve sistema reduzido: KRED * URED = FRED
      CALL GAUSS(NFREE, MAX2N, KRED, FRED, URED)

C----- Reconstroi vetor completo U
      DO I=1,NDOF
         U(I) = 0.0D0
      END DO
      DO II=1,NFREE
         I = FREE_GDLS(II)
         U(I) = URED(II)
      END DO
      DO I=1,NBC
         G1 = GDL(BCNODE(I), BCDOF(I))
         U(G1) = BCVAL(I)
      END DO

C----- Reacoes: R = K_original * U - F_original
      DO I=1,NDOF
         R(I) = 0.0D0
         DO J=1,NDOF
            R(I) = R(I) + KMAT_ORIG(I,J)*U(J)
         END DO
         R(I) = R(I) - F_ORIG(I)
      END DO

      RETURN
      END

C=============================================================
      SUBROUTINE ASSEMBLE4(KMAT, LDA, KE, G1,G2,G3,G4)
      IMPLICIT NONE
      INTEGER LDA, G1,G2,G3,G4
      DOUBLE PRECISION KMAT(LDA, *), KE(4,4)
      INTEGER I,J
      INTEGER G(4)

      G(1)=G1
      G(2)=G2
      G(3)=G3
      G(4)=G4

      DO I=1,4
         DO J=1,4
            KMAT(G(I),G(J)) = KMAT(G(I),G(J)) + KE(I,J)
         END DO
      END DO

      RETURN
      END

C=============================================================
      SUBROUTINE GAUSS(N, LDA, A, B, X)
      IMPLICIT NONE
      INTEGER N, LDA, I, J, K
      DOUBLE PRECISION A(LDA, N), B(N), X(N)
C  Resolve A*X=B por eliminacao de Gauss (sem pivotamento)
      DOUBLE PRECISION PIV, F
      DOUBLE PRECISION SUM

C----- Eliminacao
      DO K=1,N-1
         PIV = A(K,K)
         IF (DABS(PIV).LT.1.0D-30) THEN
          WRITE(*,*) 'Pivo ~ 0 em K=',K,'; sistema pode estar singular.'
          STOP
         END IF
         DO I=K+1,N
            F = A(I,K)/PIV
            A(I,K) = 0.0D0
            DO J=K+1,N
               A(I,J) = A(I,J) - F*A(K,J)
            END DO
            B(I) = B(I) - F*B(K)
         END DO
      END DO

C----- Retrosubstituicao
      X(N) = B(N)/A(N,N)
      DO I=N-1,1,-1
         SUM = 0.0D0
         DO J=I+1,N
            SUM = SUM + A(I,J)*X(J)
         END DO
         X(I) = (B(I)-SUM)/A(I,I)
      END DO

      RETURN
      END