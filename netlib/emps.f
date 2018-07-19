      PROGRAM EMPS
C
C  Program (in the Fortran 77 Subset language) to expand compressed
C  MPS-format linear-programming problems...
C
C  You may have to adjust the following assignments to STDIN, STDOUT,
C  STDERR, the Fortran unit numbers for the compressed input, expanded
C  output, and error messages (none if all goes well).  You may also
C  have to add OPEN statements for these files.
C
C  The variables X and Y declared INTEGER*4 in subroutine EXFORM must
C  carry the equivalent of at least 31 binary bits of precision.  The
C  variables X and X2 declared INTEGER*4 in subroutine CHKCHR must
C  carry at least 17 bits of precision (16 if regarded as unsigned).
C  Other integer variables can have less precision (e.g. the equivalent
C  of 16 binary bits).
C
C  Comments at the end of this program give sample input data and the
C  output it should produce.
C
C  The DATA statement a few lines below and the sample input data at
C  the end are case-sensitive.
C
C  If your Fortran I/O library does not honor carriage controls, then
C  you will find an extra blank at the beginning of each line of MPS
C  output written by this program.  To get rid of these extra blanks,
C  globally change '(1X,' into '(' in the source code.
C
      COMMON /CHKSUM/ NCS, NLINE
      INTEGER NCS, NLINE
      COMMON /ROWNAM/ RN
      CHARACTER*8 RN(8191)
      COMMON /COLNAM/ CN
      CHARACTER*8 CN(8191)
      COMMON /SSPARS/ SP
      CHARACTER*15 SP(4368)
      COMMON /STDIO/ STDERR, STOUT
      INTEGER STDERR, STOUT
      COMMON /TRNSLT/ TR0, TR
      INTEGER TR0, TR(255)
      COMMON /INVTRN/ INVTR
      CHARACTER*92 INVTR, TRINIT
      CHARACTER*1 TRI(128)
      EQUIVALENCE (TRINIT, TRI(1))
      INTEGER I, STDIN, STDOUT
      DATA TRINIT/'!"#$%&''()*+,-./0123456789;<=>?@ABCDEFGHIJKLMNOPQRSTU
     1VWXYZ[]^_`abcdefghijklmnopqrstuvwxyz{|}~'/
C
      DATA STDIN/5/, STDOUT/6/
C
C  ***  BODY  ***
C
      INVTR = TRINIT
      NLINE = 0
      TR0 = 92
      DO 10 I = 1, 255
 10      TR(I) = 92
      DO 20 I = 1, 92
 20      TR(ICHAR(TRI(I))) = I - 1
      STDERR = 0
      CALL EXPAND(STDIN,STDOUT)
      END
      SUBROUTINE EXPAND(STDIN,STDOUT)
      INTEGER STDIN, STDOUT
C
      COMMON /STDIO/ STDERR, STOUT
      INTEGER STDERR, STOUT
      COMMON /COLNAM/ CN
      CHARACTER*8 CN(8191)
      COMMON /ROWNAM/ RN
      CHARACTER*8 RN(8191)
      COMMON /SSPARS/ SP
      CHARACTER*15 SP(4368)
      COMMON /CHKSUM/ NCS, NLINE
      INTEGER NCS, NLINE
      COMMON /CHEKBF/ CHKBUF
      CHARACTER*72 CHKBUF
C
      EXTERNAL TOOBIG, EXFORM
      LOGICAL TOOBIG
C
      CHARACTER*80 LINE
      CHARACTER*72 LINE72
      CHARACTER*71 LINE2
      CHARACTER*1 L1(80), ROWTYP
      CHARACTER*4 LINE4
      CHARACTER*8 RNAME
      EQUIVALENCE (L1(1), LINE, LINE4, ROWTYP, LINE72), (L1(2), RNAME)
      EQUIVALENCE (L1(2), LINE2)
      INTEGER BDNZ, COLMX, I, J, K, NBD, NCOL, NRAN, NRHS, NROW,
     1        NS, NZ, RANZ, RHSNZ
      LOGICAL CANEND, OVRSIZ
C
C  ***  CNMAX, RNMAX, and SPMAX are the dimensions of CN, RN, and SP
C  ***  respectively.  They should be PARAMETERs, but the Fortran 77
C  ***  subset does not include PARAMETER statements.
C
      INTEGER CNMAX, RNMAX, SPMAX
      DATA CNMAX/8191/, RNMAX/8191/, SPMAX/4368/
C
C  ***  BODY  ***
C
      STOUT = STDOUT
      LINE = ' '
      CHKBUF = ' '
      NCS = 1
      CANEND = .FALSE.
C
C   *** NAME LINE ***
C
 10   READ(STDIN,'(A)',END=120) LINE
 20   NLINE = NLINE + 1
      IF (LINE4 .NE. 'NAME') GO TO 10
      CANEND = .FALSE.
      WRITE(STDOUT,'(1X,A)') LINE
C
C  *** PROBLEM STATISTICS ***
C
      READ(STDIN,'(8I9/3I9)',END=120) NROW, NCOL, COLMX, NZ, NRHS,
     1      RHSNZ, NRAN, RANZ, NBD, BDNZ, NS
      OVRSIZ = TOOBIG(NROW, RNMAX, 'RNMAX') .OR.
     1         TOOBIG(NCOL, CNMAX, 'CNMAX') .OR.
     2         TOOBIG(NS,   SPMAX, 'SPMAX')
      IF (OVRSIZ) STOP
      NLINE = NLINE + 2
C
C  ***  READ AND EXPAND NUMBER TABLE ***
C
      I = 0
      J = 73
 30   IF (I .GE. NS) GO TO 40
      I = I + 1
      IF (L1(J) .EQ. ' ') CALL RDLINE(J, LINE72, STDIN)
      CALL EXFORM(LINE, J, SP(I))
      GO TO 30
C
C  ***  READ, PRINT ROW NAMES ***
C
 40   WRITE(STDOUT,'(1X,''ROWS'')')
      DO 50 I = 1, NROW
         CALL RDLINE(J, LINE72, STDIN)
         RN(I) = RNAME
         WRITE(STDOUT, '(1X,1X,A1,2X,A8)') ROWTYP, RNAME
 50      CONTINUE
C
C  *** READ, PRINT COLUMNS ***
C
      WRITE(STDOUT,'(1X,''COLUMNS'')')
      CALL COLOUT(NZ, RN, STDIN, STDOUT, 1)
C
C  *** RIGHT-HAND SIDES ***
C
      WRITE(STDOUT,'(1X,''RHS'')')
      CALL COLOUT(RHSNZ, RN, STDIN, STDOUT, 2)
C
C  *** RANGES ***
C
      IF (NRAN .LE. 0) GO TO 60
      WRITE(STDOUT,'(1X,''RANGES'')')
      CALL COLOUT(RANZ, RN, STDIN, STDOUT, 3)
C
C  *** BOUNDS ***
C
 60   IF (NBD .LE. 0) GO TO 70
      WRITE(STDOUT,'(1X,''BOUNDS'')')
      CALL COLOUT(BDNZ, CN, STDIN, STDOUT, 4)
C
 70   IF (NCS .EQ. 1) GO TO 110
 80   READ(STDIN,'(A72)',END=120) LINE72
      NLINE = NLINE + 1
      IF (CHKBUF .NE. LINE72) THEN
         IF (L1(1) .EQ. ':') THEN
            DO 90 K = 72, 1, -1
               IF (L1(K) .NE. ' ') GO TO 100
 90            CONTINUE
            WRITE(STDERR,*) 'LINE ', NLINE,' IS EMPTY!'
            STOP
 100        CALL CHKCHR(K, L1)
            WRITE(STDOUT, '(1X,A)') LINE2
            GO TO 80
            ENDIF
         CALL BADCHK(LINE72,STDERR)
         ENDIF
 110  WRITE(STDOUT,'(1X,''ENDATA'')')
      READ(STDIN,'(A)',END=999) LINE
      CANEND = .TRUE.
      CHKBUF = ' '
      NCS = 1
      GO TO 20
 120  IF (CANEND) GO TO 999
      WRITE(STDERR,*) 'Premature end of file after line ', NLINE
      STOP
 999  RETURN
      END
      SUBROUTINE COLOUT(NZ, RN, STDIN, STDOUT, WHAT)
      INTEGER NZ, STDIN, STDOUT, WHAT
      CHARACTER*8 RN(*)
C
      COMMON /COLNAM/ CN
      CHARACTER*8 CN(8191)
      COMMON /STDIO/ STDERR, STOUT
      INTEGER STDERR, STOUT
C
      EXTERNAL EXFORM, EXINDX
      INTEGER EXINDX
C
      CHARACTER*80 LINE
      CHARACTER*72 LINE72
      CHARACTER*1 CC1(8), L1(80)
      CHARACTER*8 CURCOL, ROWNM(2)
      CHARACTER*15 RCOEF(2)
      EQUIVALENCE (L1(1), LINE, LINE72),  (CC1(1), CURCOL)
      INTEGER I, J, K, L, M, N
      CHARACTER*2 BT(6)
C
      DATA BT(1)/'UP'/, BT(2)/'LO'/, BT(3)/'FX'/, BT(4)/'FR'/,
     1     BT(5)/'MI'/, BT(6)/'PL'/
C
C  ***  BODY  ***
C
      J = 73
      I = 0
      K = 1
      M = 0
      CURCOL = ' '
      LINE = ' '
 10   IF (I .GE. NZ) GO TO 120
      I = I + 1
      IF (L1(J) .NE. ' ') GO TO 30
 20      CALL RDLINE(J, LINE72, STDIN)
 30   N = EXINDX(LINE, J)
      IF (N .NE. 0) GO TO 50
      IF (K .EQ. 2) WRITE(STDOUT,90) CURCOL, ROWNM(1), RCOEF(1)
      K = 1
      CURCOL = ' '
      DO 40 L = 1, 8
         CC1(L) = L1(J)
         J = J + 1
 40      CONTINUE
      M = M + 1
      IF (WHAT .EQ. 1) CN(M) = CURCOL
      GO TO 20
C
 50   IF (WHAT .LT. 4) GO TO 70
      IF (N .LT. 7) GO TO 60
         WRITE(STDERR,*) 'BAD BOUND TYPE INDEX =', N
         STOP
 60   IF (L1(J) .EQ. ' ') CALL RDLINE(J, LINE72, STDIN)
      ROWNM(1) = RN(EXINDX(LINE,J))
      IF (N .LT. 4) GO TO 80
         WRITE(STDOUT,110) BT(N), CURCOL, ROWNM(1)
         GO TO 10
 70   ROWNM(K) = RN(N)
 80   IF (L1(J) .EQ. ' ') CALL RDLINE(J, LINE72, STDIN)
      CALL EXFORM(LINE, J, RCOEF(K))
      IF (WHAT .GT. 3) GO TO 100
      K = K + 1
      IF (K .EQ. 2) GO TO 10
      WRITE(STDOUT,90) CURCOL, (ROWNM(K), RCOEF(K), K = 1, 2)
 90   FORMAT(1X,4X,A8,2X,A8,2X,A15,A8,2X,A15)
      K = 1
      GO TO 10
C
 100  WRITE(STDOUT,110) BT(N), CURCOL, ROWNM(1), RCOEF(1)
 110  FORMAT(1X,1X,A2,1X,A8,2X,A8,2X,A15)
      GO TO 10
 120  IF (K .EQ. 2) WRITE(STDOUT,90) CURCOL, ROWNM(1), RCOEF(1)
 999   RETURN
      END
      SUBROUTINE CHKCHR(K, L1)
      INTEGER K
      CHARACTER*1 L1(72)
      INTEGER*4 X, X2
      INTEGER I
      COMMON /CHEKBF/ CHKBUF
      CHARACTER*72 CHKBUF
      COMMON /CHKSUM/ NCS, NLINE
      INTEGER NCS, NLINE
      COMMON /TRNSLT/ TR0, TR
      INTEGER TR0, TR(255)
      COMMON /INVTRN/ INVTR
      CHARACTER*92 INVTR
      CHARACTER*1 C1(72), INTR(92)
      EQUIVALENCE (INVTR, INTR(1)), (C1(1), CHKBUF)
      X = 0
      DO 10 I = 1, K
         X2 = X / 2
         IF (2*X2 .LT. X) X2 = X2 + 16384
         X = X2 + TR(ICHAR(L1(I)))
 10      CONTINUE
      NCS = NCS + 1
      C1(NCS) = INTR(MOD(X,92)+1)
      END
      SUBROUTINE RDLINE(J, OLINE, STDIN)
      INTEGER J, STDIN
      CHARACTER*72 OLINE
C
      COMMON /STDIO/ STDERR, STOUT
      INTEGER STDERR, STOUT
      COMMON /CHKSUM/ NCS, NLINE
      INTEGER NCS, NLINE
      COMMON /CHEKBF/ CHKBUF
      CHARACTER*72 CHKBUF
C
      INTEGER K
      CHARACTER*72 CHKLIN, LINE
      CHARACTER*71 LINE2
      CHARACTER*1 C1(72), CL1(72), L1(72)
      EQUIVALENCE (LINE, L1(1)), (C1(1), CHKBUF), (CL1(1), CHKLIN)
      EQUIVALENCE (L1(2), LINE2)
C
C  ***  BODY  ***
C
 10   READ(STDIN,'(A72)',END=40) LINE
      OLINE = LINE
      J = 1
      NLINE = NLINE + 1
      DO 20 K = 72, 1, -1
         IF (L1(K) .NE. ' ') GO TO 30
 20      CONTINUE
      WRITE(STDERR,*) 'LINE ', NLINE,' IS EMPTY!'
      STOP
 30   CALL CHKCHR(K, L1)
      IF (NCS .LT. 72) GO TO 50
      READ(STDIN,'(A72)',END=40) CHKLIN
      NLINE = NLINE + 1
      IF (CHKBUF .NE. CHKLIN) CALL BADCHK(CHKLIN,STDERR)
      CHKBUF = ' '
      NCS = 1
      GO TO 50
 40   WRITE(STDERR,*) 'Premature end of file after line ', NLINE
      STOP
 50   IF (L1(1) .EQ. ':') THEN
         WRITE(STOUT, '(1X,A)') LINE2
         GO TO 10
         ENDIF
      END
      SUBROUTINE BADCHK(BLINE,STDERR)
      CHARACTER*72 BLINE
      INTEGER STDERR
C
      COMMON /CHKSUM/ NCS, NLINE
      INTEGER NCS, NLINE
      COMMON /CHEKBF/ CHKBUF
      CHARACTER*72 CHKBUF
C
      INTEGER I
      CHARACTER*72 LINE
      CHARACTER*1 C1(72), CL1(72)
      EQUIVALENCE (C1(1), CHKBUF), (CL1(1), LINE)
C
      LINE = BLINE
      WRITE(STDERR,10) CHKBUF, LINE
 10   FORMAT(' Check sum error: expected'/1X,A/' but got'/1X,A)
      IF (CL1(1) .EQ. ' ') GO TO 20
         WRITE(STDERR,*) 'Line ', NLINE, ' = bad check sum line'
         STOP
 20      I = 1
 30      I = I + 1
         IF (CL1(I) .EQ. C1(I)) GO TO 30
         WRITE(STDERR,*) 'Bad check sum for line ', I + NLINE - NCS - 1
         STOP
      END
      INTEGER FUNCTION EXINDX(L1, J)
      CHARACTER*1 L1(72)
      INTEGER J
C
      COMMON /STDIO/ STDERR, STOUT
      INTEGER STDERR, STOUT
      COMMON /TRNSLT/ TR0, TR
      INTEGER TR0, TR(255)
      COMMON /CHKSUM/ NCS, NLINE
      INTEGER NCS, NLINE
C
      INTEGER I, K, X
C
C  ***  BODY  ***
C
      K = TR(ICHAR(L1(J)))
      IF (K .LT. 46) GO TO 10
      WRITE(STDERR,*) 'Error on line ', NLINE,
     1                ' -- EXINDX not given index; arg ='
      WRITE(STDERR,'(1X,72A1)') (L1(I), I = J, 72)
      STOP
 10   J = J + 1
      IF (K .LT. 23) GO TO 20
         EXINDX = K - 23
         GO TO 999
 20   X = K
 30   K = TR(ICHAR(L1(J)))
      J = J + 1
      X = X*46 + K
      IF (K .LT. 46) GO TO 30
      X = X - 46
      EXINDX = X
 999  RETURN
      END
      SUBROUTINE STORE(A, B)
      CHARACTER*15 A, B
      A = B
      END
      SUBROUTINE EXFORM(L1, J, RC)
      CHARACTER*1 L1(72), RC(15)
      INTEGER J
C
      COMMON /SSPARS/ SP
      CHARACTER*15 SP(4368)
      COMMON /TRNSLT/ TR0, TR
      INTEGER TR0, TR(255)
C
      EXTERNAL EXINDX
      INTEGER EXINDX
C
      INTEGER EX, I, K, L, ND, NELIM
C  ***  X NEEDS AT LEAST 31 BITS OF PRECISION  ***
      INTEGER*4 X, Y
      CHARACTER*1 D(15), DIGITS(10)
      CHARACTER*10 DIG10
      EQUIVALENCE(DIG10, DIGITS(1))
      DATA DIG10/'0123456789'/
C
C  ***  BODY  ***
C
      K = TR(ICHAR(L1(J)))
      IF (K .GE. 46) GO TO 10
C
C  ***  SUPER-SPARSE ***
C
         CALL STORE(RC, SP(EXINDX(L1,J)))
         GO TO 999
 10   J = J + 1
      K = K - 46
      IF (K .LT. 23) GO TO 20
         K = K - 23
         RC(1) = '-'
         L = 2
         NELIM = 11
         GO TO 30
 20   L = 1
      NELIM = 12
 30   IF (K .LT. 11) GO TO 90
C
C  ***  INTEGER FLOATING-POINT ***
C
      K = K - 11
      D(1) = '.'
      I = 1
      IF (K .LT. 6) GO TO 40
         X = K - 6
         GO TO 60
 40   X = K
 50      K = TR(ICHAR(L1(J)))
         J = J + 1
         X = X*46 + K
         IF (K .LT. 46) GO TO 50
      X = X - 46
 60   IF (X .NE. 0) GO TO 70
         D(2) = '0'
         I = 2
         GO TO 80
 70   I = I + 1
      D(I) = DIGITS(1 + MOD(X,10))
      X = X / 10
      IF (X .NE. 0) GO TO 70
 80      RC(L) = D(I)
         L = L + 1
         I = I - 1
         IF (I .GT. 0) GO TO 80
      GO TO 370
C
C  ***  NON-INTEGER FLOATING-POINT ***
C
 90   EX = TR(ICHAR(L1(J))) - 50
      X = TR(ICHAR(L1(J+1)))
      J = J + 2
      Y = 0
      I = 0
 100  IF (K .LE. 0) GO TO 130
         K = K - 1
         IF (X .LT. 100000000) GO TO 110
            Y = X
            X = TR(ICHAR(L1(J)))
            GO TO 120
 110     X = X*92 + TR(ICHAR(L1(J)))
 120     J = J + 1
         GO TO 100
 130  IF (Y .EQ. 0) GO TO 150
 140     I = I + 1
         D(I) = DIGITS(1 + MOD(X,10))
         X = X / 10
         IF (X .GT. 1) GO TO 140
      X = Y
      GO TO 160
 150  IF (X .NE. 0) GO TO 160
            D(1) = '0'
            I = 1
            GO TO 170
 160     I = I + 1
         D(I) = DIGITS(1 + MOD(X,10))
         IF (X .LT. 10) GO TO 170
         X = X / 10
         GO TO 160
 170  ND = I + EX
      IF (EX .LE. 0) GO TO 210
      IF (ND .LT. NELIM) GO TO 180
      IF (EX .GE. 3) GO TO 280
 180     RC(L) = D(I)
         L = L + 1
         I = I - 1
         IF (I .GE. 1) GO TO 180
C
 190     RC(L) = '0'
         L = L + 1
         IF (EX .LE. 1) GO TO 200
         EX = EX - 1
         GO TO 190
 200  RC(L) = '.'
      L = L + 1
      GO TO 370
 210  IF (ND .LT. 0) GO TO 250
 220     ND = ND - 1
         IF (ND .LT. 0) GO TO 230
         RC(L) = D(I)
         L = L + 1
         I = I - 1
         GO TO 220
 230  RC(L) = '.'
      L = L + 1
 240     IF (I .LE. 0) GO TO 370
         RC(L) = D(I)
         L = L + 1
         I = I - 1
         GO TO 240
 250  IF (EX .LE. -NELIM) GO TO 280
         RC(L) = '.'
         L = L + 1
 260     IF (ND .GE. 0) GO TO 270
            ND = ND + 1
            RC(L) = '0'
            L = L + 1
            GO TO 260
 270     IF (I .LE. 0) GO TO 370
            RC(L) = D(I)
            L = L + 1
            I = I - 1
            GO TO 270
 280  EX = EX + I - 1
      IF (-10 .NE. EX) GO TO 290
         EX = -9
         GO TO 320
 290  IF (EX .LE. 9) GO TO 310
      IF (EX - I .GT. 8) GO TO 310
 300     RC(L) = D(I)
         L = L + 1
         I = I - 1
         EX = EX - 1
         IF (EX .GT. 9) GO TO 300
 310  RC(L) = D(I)
      I = I - 1
      L = L + 1
 320  RC(L) = '.'
      L = L + 1
 330     IF (I .LE. 0) GO TO 340
            RC(L) = D(I)
            L = L + 1
            I = I - 1
            GO TO 330
 340     RC(L) = 'E'
         L = L + 1
         IF (EX .GE. 0) GO TO 350
            RC(L) = '-'
            L = L + 1
            EX = -EX
 350     IF (EX .LE. 0) GO TO 360
            I = I + 1
            D(I) = DIGITS(1 + MOD(EX,10))
            EX = EX / 10
            GO TO 350
 360     IF (I .LE. 0) GO TO 370
            RC(L) = D(I)
            L = L + 1
            I = I - 1
            GO TO 360
 370  CONTINUE
C  ***  COMMENT OUT THE NEXT 10 LINES TO GET LEFT-JUSTIFIED NUMBERS ***
      IF (L .GT. 12) GO TO 400
         I = 12
 380        L = L - 1
            RC(I) = RC(L)
            I = I - 1
            IF (L .GT. 1) GO TO 380
 390     RC(I) = ' '
         I = I - 1
         IF (I .GE. 1) GO TO 390
      L = 13
 400  IF (L .GT. 15) GO TO 999
         RC(L) = ' '
         L = L + 1
         GO TO 400
 999  RETURN
      END
      LOGICAL FUNCTION TOOBIG(N, NMAX, WHAT)
      INTEGER N, NMAX
      CHARACTER*5 WHAT
C
      COMMON /STDIO/ STDERR, STOUT
      INTEGER STDERR, STOUT
C
C  ***  BODY  ***
C
      TOOBIG = .FALSE.
      IF (N .LE. NMAX) GO TO 999
      TOOBIG = .TRUE.
      WRITE(STDERR,'('' Increase '',A5,'' from '',I10,'' to at least '',
     1I10,'' and try again!'')') WHAT, NMAX, N
 999  RETURN
C
C Test data:  if you remove the C in column 1 of the following
C 30 lines of compressed data (so that NAME begins in column 1 and
C the final 6 of the line after NAME is in column 72)....
C
CNAME          EMPSTEST
C        7        8        2       16        2        5        2        6
C        1        8        5
CSO"$KYiR#TssQR#TiQ,|
CEVLRES
CNOBJEC
CERAI72
CGDEP73
CLDEP72
CETRS72
CGINV72
C8RVAD72
C<c;c8RVAD73
C<QQ,|;c8DEPN72
C>>=c8DEPN73
C<>?c8INVT72
C@z?z8WK1T78
C@z?z8WK2T78
C@{?z8WK3T78
C@|?z
C8RHS1
C<iQ.H?QR'0@{8RHS2
C@}?g
C8RAN1
C9=<<>;8RAN2
C9=<<>;
C8BNDS1
C=9?;9?9>@;ARO3_x9<RO2{"<=9<>9
C .BZO1nPy"=C75]/DMK'>CzASq
C
C and feed the result to the above program, then the following (minus
C the C's in column 1) is the output you should get...
CNAME          EMPSTEST
CROWS
C E  VLRES
C N  OBJEC
C E  RAI72
C G  DEP73
C L  DEP72
C E  TRS72
C G  INV72
CCOLUMNS
C    RVAD72    RAI72               1.   OBJEC               1.
C    RVAD73    RAI72            1.101   OBJEC               1.
C    DEPN72    DEP72           -1.101   DEP73               1.
C    DEPN73    RAI72           -1.101   TRS72               1.
C    INVT72    INV72              -1.   TRS72              -1.
C    WK1T78    INV72              -1.   TRS72              -1.
C    WK2T78    INV72              -2.   TRS72              -1.
C    WK3T78    INV72              -3.   TRS72              -1.
CRHS
C    RHS1      RAI72           -1.234   TRS72             5.67
C    RHS1      INV72              -2.
C    RHS2      INV72              -4.   TRS72               5.
CRANGES
C    RAN1      VLRES             2.34   RAI72             -34.
C    RAN1      DEP72            -2.34
C    RAN2      VLRES             2.34   RAI72             -34.
C    RAN2      DEP72            -2.34
CBOUNDS
C FR BNDS1     RVAD72
C PL BNDS1     RVAD73
C UP BNDS1     WK1T78         8.07907
C MI BNDS1     WK2T78
C LO BNDS1     WK3T78         1.57957
C UP BNDS1     DEPN72         1.51985
C FX BNDS1     DEPN73         8.07907
C FX BNDS1     INVT72         8.07907
CENDATA
      END
