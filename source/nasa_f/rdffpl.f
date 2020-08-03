C$Procedure   RDFFPL ( read triangular plate model flatfile )
 
      SUBROUTINE RDFFPL ( INFILE, PLTTYP, NV, VRTCES, NP, PLATES )
      IMPLICIT NONE
 
C$ Abstract
C
C     Read a triangular plate model's flatfile's data file.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     MKPLAT User's Guide
C     SURFACE
C
C$ Keywords
C
C     PLATE
C     FLATFILE
C
C$ Declarations
  
      INCLUDE               'mkdsk02.inc'
      INCLUDE               'mkdsk.inc'

      CHARACTER*(*)         INFILE
      INTEGER               PLTTYP 
      INTEGER               NV
      DOUBLE PRECISION      VRTCES ( 3, * )
      INTEGER               NP
      INTEGER               PLATES ( 3, * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFILE     I   Name of plate data file.
C     PLTTYP     I   Integer code identifying type of data file.
C     NV         O   Number of vertices in model.
C     VRTCES     O   3 x NV array of vertices.
C     NP         O   Number of plates in model.
C     PLATES     O   3 x NP array of plates.
C
C
C$ Detailed_Input
C
C     INFILE     is the full pathname of the plate model flatfile.
C
C     PLTTYP     the integer code identifying the type of data
C                to read from INFILE. Two values are allowed:
C
C                    1    A standard plate-vertex data file
C                         containing vertex coordinates and
C                         the plate-vertex mappings.
C
C                    2    A Gaskel shape data file containing
C                         plate ordered vertex data.
C
C                    3    A vertex-facet table data file
C                         containing vertex coordinates and
C                         the facet listing.
C
C                    4    Rosetta/Osiris ".ver" format file
C                         containing vertex coordinates
C                         and a plate-vertex mapping.
C
C$ Detailed_Output
C
C     NV         is the number of vertices in the plate model.
C
C
C     VRTCES     is an array of the NV vertices given in a body-fixed
C                frame of reference. Elements 
C
C                   VRTCES(J,I), J = 1 ... 3
C                      
C                are, respectively, the X, Y, and Z coordinates of the
C                Ith vertex.
C
C
C     PLATES     is an array of the NP plates. Elements 
C
C                   PLATES(J,I), J = 1 ... 3
C                      
C                are, respectively, the indices of the vertices of the
C                Ith plate. The vertex indices range from 1 to NV. 
C
C                The order of the vertices give the direction of the
C                plate's outward normal vector: the vector is parallel
C                to the cross product of the plate edge connecting
C                vertex 1 to vertex 2 and the plate edge connecting
C                vertex 2 to vertex 3, in that order.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the number of vertices (NV) exceeds the value of MAXVRT,
C        SPICE(TOOMANYVERTICES) signals.
C
C     2) If the number of plates (NP) exceeds the value of MAXPLT,
C        SPICE(TOOMANYPLATES) signals.
C
C     3) If the NPARSD routine fails while parsing a double string, 
C        SPICE(BADDOUBLEPRECISION) signals.
C
C     4) If the NPARSI routine fails while parsing an integer string, 
C        SPICE(BADINTEGER) signals.
C
C     5) If a data line lacks the expected format (data type 3),
C        SPICE(BADDATALINE) signals.
C
C     6) If the PLTTYP variable has an uncoded value,
C        SPICE(BADDATATYPE) signals.
C
C     7) If the Gaskell ICQ parameter "Q" is too large, the error
C        SPICE(QPARAMOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     Input File Format
C     -----------------
C
C     INFILE is assumed to be structured as follows (values
C            are space delimited within each input record).
C
C     -Plate/vertex data file (PLATE_TYPE 1):
C
C        Line  1: NV
C        NV = the number of vertices (INTEGER)
C
C        Lines 2 to (NV+1): VID X Y Z
C        VID         = the vertex ID (INTEGER)
C        X, Y, and Z = the vertex coordinates (DOUBLE PRECISION)
C
C        Line  (NV+2): NP
C        NP = number of triangular plates (INTEGER)
C
C        Lines (NV+3) to (NV+3+NP+1): PID V1 V2 V3
C        PID            = the plate ID (INTEGER)
C        V1, V2, and V3 = the plate vertex ID's (INTEGER)
C
C     -Shape data file (PLATE_TYPE 2):
C
C        Line 1: Q
C        Q = number of vertex coordinates per side of the 
C            cube face - not counting the origin. A
C            face has (Q+1)^2 vertex coordinates.
C
C        Lines 2 to (Q+1): X Y Z
C        X, Y, and Z = the vertex coordinates (DOUBLE PRECISION)
C
C     -Vertex/Facet table data file (PLATE_TYPE 3):
C
C        NV = the number of vertices (INTEGER)
C
C        Lines 1 to NV: 'V' X Y Z
C        'V'         = character flag indicating vertex data
C        X, Y, and Z = the vertex coordinates (DOUBLE PRECISION)
C
C        NP = number of triangular plates (INTEGER)
C
C        Lines (NV+1) to (NV+1+NP): 'F' V1 V2 V3
C        'F'            = character flag indicating facet 
C                         (plate-vertex) data
C        V1, V2, and V3 = the plate vertex ID's (INTEGER)
C
C$ Examples
C
C  C     The following include file defines MAXVRT and
C  C     MAXPLT, the maximum number of vertices and plates
C  C     used in the plate model software.
C
C        INCLUDE               'mkdsk02.inc'
C
C        INTEGER               NV
C        INTEGER               NP
C        INTEGER               PLATES ( 3, MAXPLT )
C        INTEGER               PLTTYP
C
C        DOUBLE PRECISION      VRTCES ( 3, MAXVRT )
C
C        CHARACTER*(80)        INFILE
C
C        CALL PROMPT ( 'Enter name of plate data file : ', INFILE )
C        CALL PROMPT ( 'Enter type of plate data      : ', PLTTYP )
C
C        CALL RDFFPL ( INFILE,  PLTTYP, NV, VERTCES, NP, PLATES )
C
C$ Restrictions
C
C     It is the user's responsibility to properly dimension
C     arrays in the calling routine large enough to accept data
C     from the given file.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     J.A. Bytof      (JPL)
C
C$ Version
C
C-    MKDSK Version 5.1.0, 21-MAR-2017 (NJB)
C
C        Now references non-portable include file 
C
C           mkdsk02.inc
C
C        instead of
C
C           dsk02.inc
C
C        This enables use of smaller maximum vertex and plate
C        counts on platforms where the default values consume
C        excessive memory.
C
C        Added checks for excessive vertex and plate counts.
C   
C        Added FAILED checks following RDFFDI, RDNBL, and
C        PRSINT calls.
C
C-    MKDSK Version 5.0.0, 25-APR-2016 (NJB)
C     
C        Bug fix: added check-out and return to input format 2 
C        branch, in the case where a Q value out of range is found.
C        
C      
C
C        07-APR-2015 (NJB)
C
C           Increased string length used to parse numeric tokens.
C           Updated message indicating completion of input file read.
C           Cleaned up debugging code.
C
C        05-AUG-2014 (NJB)
C
C           Argument list change: arguments PID and VID have been
C           removed; arguments X, Y, Z and V1, V2, V3 have been
C           replaced by arrays VRTCES and PLATES respectively.
C
C           Local array VEC is no longer used.
C     
C-    MKDSK Version 4.0.0, 08-JUN-2010 (NJB)
C
C        Added capability of reading Rosetta Osiris ".ver"
C        format file.
C
C-    MKDSK Version 3.1.0, 04-MAY-2010 (NJB)
C
C        Changed INCLUDE file from platmax.inc to dsk02.inc. 
C        Added INCLUDE statement referencing mkdsk.inc to 
C        declare parameter MAXQ.
C
C-    MKDSK Version 3.0.1, 08-OCT-2009 (NJB)
C
C        Re-ordered header sections.
C
C-    MKDSK Version 3.0.0, 25-OCT-2004 (EDW)
C
C        Added capability to process Bob Gaskell shape files.
C        Added MAXQ parameter to 'pltmax.inc'.
C
C-    MKDSK Version 2.0.0, 22-OCT-1998 (JAB)
C
C-    MKDSK Version 1.0.0, 09-APR-1997 (JAB)
C
C-&
 
C$ Index_Entries
C
C     read triangular plate model flatfile
C
C-&
 
C
C     SPICELIB functions
C 
      DOUBLE PRECISION      VNORM

      LOGICAL               EQSTR
      LOGICAL               EVEN
      LOGICAL               FAILED
      LOGICAL               RETURN
 
C
C     Local parameters
C
      INTEGER               MAXWDS
      PARAMETER           ( MAXWDS = 3 )

      INTEGER               SHORT
      PARAMETER           ( SHORT = 80 )

      INTEGER               TOKLEN
      PARAMETER           ( TOKLEN = 40 )
C
C     Local variables
C 
      CHARACTER*(TOKLEN)    ARC     ( 6 )
      CHARACTER*(160)       ERROR
      CHARACTER*(10)        FORMAT
      CHARACTER*(LNSIZE)    ILINE     
      CHARACTER*(SHORT)     TOKENS ( MAXWDS )

      DOUBLE PRECISION      ARD    ( 3 )
      DOUBLE PRECISION      VEC    ( 3 )
      DOUBLE PRECISION      W1     ( 3 )
      DOUBLE PRECISION      W2     ( 3 )

      INTEGER               ARI    ( 4 )
      INTEGER               FACE
      INTEGER               I
      INTEGER               IX1
      INTEGER               IX2
      INTEGER               J
      INTEGER               N      ( 0:MAXQ, 0:MAXQ, 6)
      INTEGER               N0
      INTEGER               NC
      INTEGER               ND
      INTEGER               NI
      INTEGER               NREC
      INTEGER               NTOK
      INTEGER               PTR 
      INTEGER               Q

      LOGICAL               EOF

C
C     Saved variables
C
      SAVE

C
C     Standard SPICE error handling.
C
 
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'RDFFPL' )

C
C     Initialize record counter.
C
      NREC = 1

C
C     Read the data file corresponding to the PLTTYP ID.
C
      IF( PLTTYP .EQ. 1 ) THEN

C
C        Plate data type 1.
C
C        First line, a single integer, number of vertices.
C
         FORMAT = 'I'
  
         CALL RDFFDI ( INFILE, NREC, FORMAT, ND, NI, NC,
     .                 ARD, ARI, ARC, EOF  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'RDFFPL' )
            RETURN
         END IF

         NV = ARI(1)

C
C        Check the routine can handle the vertices.
C
         IF ( NV .GT. MAXVRT ) THEN
            CALL SETMSG ( 'Number of vertices # exceeds limit #.' )
            CALL ERRINT ( '#',  NV                                )
            CALL ERRINT ( '#',  MAXVRT                            )
            CALL SIGERR ( 'SPICE(TOOMANYVERTICES)'                )
            CALL CHKOUT ( 'RDFFPL'                                )
            RETURN
         END IF

C
C        Read data: integer double double double.
C
         FORMAT = 'I D D D'

         DO I = 1, NV
      
            NREC   = NREC + 1
         
            CALL RDFFDI ( INFILE, NREC, FORMAT, ND, NI, NC,
     .                    ARD, ARI, ARC, EOF  )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'RDFFPL' )
               RETURN
            END IF

            IF ( EOF ) THEN
               CALL EXTMSI ( 'End of file after line #.', '#', NREC )
            END IF

            VRTCES(1,I) = ARD(1)
            VRTCES(2,I) = ARD(2)
            VRTCES(3,I) = ARD(3)
   
         END DO

C
C        Plate data type 1 reads the plat-vertex information from
C        the plate file. 
C
C        Read number of plates.
C
         NREC   = NREC + 1
         FORMAT = 'I'
 
         CALL RDFFDI ( INFILE, NREC, FORMAT, ND, NI, NC,
     .                 ARD, ARI, ARC, EOF  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'RDFFPL' )
            RETURN
         END IF

         IF ( EOF ) THEN
            CALL EXTMSI ( 'End of file after line #.', '#', NREC )
         END IF
 
C
C        Check we can process 'NP' plates.
C
         NP = ARI(1)

         IF ( NP .GT. MAXPLT ) THEN
            CALL SETMSG ( 'Number of plates # exceeds limit #.' )
            CALL ERRINT ( '#', NP                               )
            CALL ERRINT ( '#', MAXPLT                           )
            CALL SIGERR ( 'SPICE(TOOMANYPLATES)'                )
            CALL CHKOUT ( 'RDFFPL'                              )
            RETURN
         END IF
 
C
C        Read each plate's ID and corresponding vertex set defining
C        the plate.
C
         FORMAT = 'I I I I'
 
         DO I = 1, NP
            
            NREC   = NREC + 1
            
            CALL RDFFDI ( INFILE, NREC, FORMAT, ND, NI, NC,
     .                    ARD, ARI, ARC, EOF  )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'RDFFPL' )
               RETURN
            END IF

            IF ( EOF ) THEN
               CALL EXTMSI ( 'End of file after line #.', '#', NREC )
            END IF

            
            PLATES(1,I) = ARI(2)
            PLATES(2,I) = ARI(3)
            PLATES(3,I) = ARI(4)

         END DO

      ELSE IF( PLTTYP .EQ. 2 ) THEN

C
C        Plate data type 2.
C
C        Read a Gaskell shape file. Shape files contain only vertex
C        data, the plate-vertex mappings implicitly known from the 
C        ordering of the data. The Gaskell model uses quadralaterals
C        plates, a plate defined by the vector set:
C
C        [ VEC( *, I, J  , FACE), VEC( *, I+1, J  , FACE), 
C          VEC( *, I, J+1, FACE), VEC( *, I+1, J+1, FACE) ]
C
C        The file lists the vertex data in terms cube faces. Six
C        faces, each with Q+1 x Q+1 vertices.
C
         FORMAT = 'I'
  
         CALL RDFFDI ( INFILE, NREC, FORMAT, ND, NI, NC,
     .                 ARD, ARI, ARC, EOF  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'RDFFPL' )
            RETURN
         END IF

         Q  = ARI(1)
         NV = 6*( (Q+1)**2)

C
C        Does the read Q exceed the maximum value.
C
         IF ( Q .GT. MAXQ ) THEN

            CALL SETMSG ( 'Shape parameter Q = #, exceeds maximum '
     .      //            'value #. This error may indicate that '
     .      //            'the input file format is something other '
     .      //            'than the Gaskell ICQ format.'              )
            CALL ERRINT ( '#', Q                                      )
            CALL ERRINT ( '#', MAXQ                                   )
            CALL SIGERR ( 'SPICE(QPARAMOUTOFRANGE)'                   )
            CALL CHKOUT ( 'RDFFPL'                                    )
            RETURN

         END IF
C
C        Check the routine can handle the vertices.
C
         IF ( NV .GT. MAXVRT ) THEN
            CALL SETMSG ( 'Number of vertices # exceeds limit #.' )
            CALL ERRINT ( '#', NV                                 )
            CALL ERRINT ( '#', MAXVRT                             )
            CALL SIGERR ( 'SPICE(TOOMANYVERTICES)'                )
            CALL CHKOUT ( 'RDFFPL'                                )
            RETURN
         END IF

C
C        Read data: double double double
C
         FORMAT = 'D D D'

         N0 = 0

         DO FACE = 1, 6

            DO J = 0, Q

               DO I = 0, Q

                  NREC = NREC + 1

C
C                 Read the vertex data. The 3-Vector
C                    _                 _
C                   |  VEC(1,I,J,FACE)  |
C                   |  VEC(2,I,J,FACE)  |
C                   |_ VEC(3,I,J,FACE) _|
C
C                 contains the vertex coordinates.
C
                  CALL RDFFDI ( INFILE, NREC, FORMAT, 
     .                          ND,     NI,   NC, 
     .                          VEC,    ARI,  ARC,
     .                          EOF )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'RDFFPL' )
                     RETURN
                  END IF


                  IF ( EOF ) THEN
                     CALL EXTMSI ( 'End of file after line #.', 
     .                             '#', NREC)
                  END IF

C
C                Store the vertex ID to the 'N' array using the same
C                indexing as 'VEC'.
C
                 N0          = N0+1
                 N(I,J,FACE) = N0
                 
C
C                Copy the vertex coordinates to the corresponding
C                plate model arrays. This operation is somewhat
C                repetitive, but clarity is most important.
C
                 VRTCES(1,N0) = VEC(1)
                 VRTCES(2,N0) = VEC(2)
                 VRTCES(3,N0) = VEC(3)

               END DO

            END DO

         END DO

C
C        Plate data type 2 calculates the plate-vertex mappings 
C        from the ordering of the vertex data.

C
C        Number of plates based on the 'Q' value.
C
         NP = 12*(Q**2)

C
C        Check we can process 'NP' plates.
C
         IF ( NP .GT. MAXPLT ) THEN
            CALL SETMSG ( 'Number of plates # exceeds limit #.' )
            CALL ERRINT ( '#', NP                               )
            CALL ERRINT ( '#', MAXPLT                           )
            CALL SIGERR ( 'SPICE(TOOMANYPLATES)'                )
            CALL CHKOUT ( 'RDFFPL'                              )
            RETURN
         END IF

C
C        Expand a Gaskell shape model vertex ordering to
C        form a plate-vertex map.
C

C
C        Cube edges. Associate later IDs with
C        the first occurrence.
C
         DO I=1,Q-1

            N(I,Q,6) = N(Q-I, Q  , 4)
            N(I,0,6) = N(I  , Q  , 2)
            N(I,0,5) = N(Q  , Q-I, 1)
            N(I,0,4) = N(Q-I, 0  , 1)
            N(I,0,3) = N(0  , I  , 1)
            N(I,0,2) = N(I  , Q  , 1)

         END DO

         DO J=1,Q-1

            N(Q,J,6) = N(J  , Q, 5)
            N(Q,J,5) = N(0  , J, 4)
            N(Q,J,4) = N(0  , J, 3)
            N(Q,J,3) = N(0  , J, 2)
            N(0,J,6) = N(Q-J, Q, 3)
            N(0,J,5) = N(Q  , J, 2)

         END DO

C
C        Cube corners. Associate later IDs with
C        the first occurrence.
C
         N(0,0,3) = N(0,0,1)
         N(Q,0,4) = N(0,0,1)
         N(0,0,2) = N(0,Q,1)
         N(Q,0,3) = N(0,Q,1)
         N(0,0,4) = N(Q,0,1)
         N(Q,0,5) = N(Q,0,1)
         N(0,0,5) = N(Q,Q,1)
         N(Q,0,2) = N(Q,Q,1)
         N(0,0,6) = N(0,Q,2)
         N(Q,Q,3) = N(0,Q,2)
         N(0,Q,5) = N(Q,Q,2)
         N(Q,0,6) = N(Q,Q,2)
         N(Q,Q,4) = N(0,Q,3)
         N(0,Q,6) = N(0,Q,3)
         N(Q,Q,5) = N(0,Q,4)
         N(Q,Q,6) = N(0,Q,4)

         N0 = 0

         DO FACE = 1, 6 

            DO I = 0, Q-1

               DO J = 0, Q-1

C
C                 As mentioned, the Gaskell shape model uses
C                 quadrilateral plates as opposed to triangular used by
C                 the NAIF plate system. We can reduce a quadrilateral
C                 into two triangles by adding a vector connecting
C                 opposite facing vertices creating two triangles
C                 (1 and 2):
C
C                               V1 
C                               /\
C                              / 1\
C                          V4 /____\ V2
C                             \  2 /
C                              \  /
C                               \/
C                               V3
C
C                 or
C                                V1
C                               /|\
C                              / | \
C                         V4  /  |2 \ V2
C                             \ 1|  /
C                              \ | /
C                               \|/
C                                V3
C
C                 We chose the connecting vector to minimize the cross
C                 product magnitude of the connected vertex vectors.
C
C                    connection = MIN( V1 x V3, V2 X V4)
C

C
C                 Connection 1.
C

                  IX1 = (FACE-1)*((Q+1)**2) + ((Q+1)* J   ) +  I    + 1
                  IX2 = (FACE-1)*((Q+1)**2) + ((Q+1)*(J+1)) + (I+1) + 1
                  
                  CALL VCRSS ( VRTCES(1,IX1), VRTCES(1,IX2), W1 )

C
C                 Connection 2.
C

                  IX1 = (FACE-1)*((Q+1)**2) + ((Q+1)* J   ) + (I+1) + 1
                  IX2 = (FACE-1)*((Q+1)**2) + ((Q+1)*(J+1)) +  I    + 1


                  CALL VCRSS ( VRTCES(1,IX1), VRTCES(1,IX2), W2 )

C
C                 Calculate the magnitudes of the cross products;
C                 branch based on the minimum value.
C
C                 Fill two entries of the plate-vertex mapping 
C                 from the divided quadralateral.
C
                  IF( VNORM(W1) .LE. VNORM(W2) ) THEN
C
C                    Apply connection 1.
C
                     N0           = N0+1
                     PLATES(1,N0) = N(I ,  J,   FACE)
                     PLATES(2,N0) = N(I+1, J+1, FACE)
                     PLATES(3,N0) = N(I+1, J,   FACE)

                     N0           = N0+1
                     PLATES(1,N0) = N(I,   J,   FACE)
                     PLATES(2,N0) = N(I,   J+1, FACE)
                     PLATES(3,N0) = N(I+1, J+1, FACE)

                  ELSE
C
C                    Apply connection 2.
C
                     N0           = N0+1
                     PLATES(1,N0) = N(I,   J,   FACE)
                     PLATES(2,N0) = N(I,   J+1, FACE)
                     PLATES(3,N0) = N(I+1, J,   FACE)

                     N0           = N0+1
                     PLATES(1,N0) = N(I+1, J,   FACE)
                     PLATES(2,N0) = N(I,   J+1, FACE)
                     PLATES(3,N0) = N(I+1, J+1, FACE)
 
                  END IF

               END DO

            END DO

         END DO




      ELSE IF( PLTTYP .EQ. 3 ) THEN
C
C        Plate data type 3.
C
C        Read data. The data format depends on the
C        prefix marker either 'v' (vertex) or 'f' (facet).
C        Since we don't know the prefix a priori,
C        extract the four data elements, then parse.
C
         NP     = 0
         NV     = 0
         FORMAT = 'C C C C'
         EOF    = .FALSE.
   
         DO WHILE ( .NOT. EOF )

            CALL RDFFDI ( INFILE, NREC, FORMAT, ND, NI, NC,
     .                    ARD, ARI, ARC, EOF  )
     
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'RDFFPL' )
               RETURN
            END IF

            IF ( .NOT. EOF ) THEN

C
C              Branch based on the type of data, 'v' or 'f'.
C
               IF(  EQSTR( ARC(1), 'v' )  ) THEN

C
C                 Vertex data. Process a line of three doubles.
C                 The algorithm includes an implicit assumption that
C                 the data order defines the vertex IDs.  If not,
C                 the output plate files will be useless.
C
                  NV = NV + 1

                  IF ( NV .GT. MAXVRT ) THEN

                     CALL SETMSG ( 'Number of vertices # exceeds '
     .               //            'limit #.'                    )
                     CALL ERRINT ( '#',  NV                      )
                     CALL ERRINT ( '#',  MAXVRT                  )
                     CALL SIGERR ( 'SPICE(TOOMANYVERTICES)'      )
                     CALL CHKOUT ( 'RDFFPL'                      )
                     RETURN

                  END IF

                  DO J = 2, 4

                     CALL NPARSD ( ARC(J), ARD(J-1), ERROR, PTR )
                  
                     IF ( PTR .NE. 0 ) THEN
                        CALL SETMSG ( 'D.P. error (#) in line #.' )
                        CALL ERRCH  ( '#', ERROR                  )
                        CALL ERRINT ( '#', NREC                   )
                        CALL SIGERR ( 'SPICE(BADDOUBLEPRECISION)' )
                        CALL CHKOUT ( 'RDFFPL'                    )
                        RETURN
                     END IF

                  END DO

   
                 VRTCES(1,NV) = ARD(1)
                 VRTCES(2,NV) = ARD(2)
                 VRTCES(3,NV) = ARD(3)


               ELSE IF ( EQSTR( ARC(1), 'f') ) THEN

C
C                 Plate-vertex data. Process a line of three integers.
C                 The algorithm includes an implicit assumption that
C                 the data order defines the plate IDs. If not,
C                 the output plate files will be useless.
C
                  NP = NP + 1


                 IF ( NP .GT. MAXPLT ) THEN

                    CALL SETMSG( 'Number of plates # exceeds limit #.' )
                    CALL ERRINT( '#', NP                               )
                    CALL ERRINT( '#', MAXPLT                           )
                    CALL SIGERR( 'SPICE(TOOMANYPLATES)'                )
                    CALL CHKOUT( 'RDFFPL'                              )
                    RETURN

                 END IF


                  DO J =2, 4

                     CALL NPARSI ( ARC(J), ARI(J-1), ERROR, PTR )
                  
                     IF ( PTR .NE. 0 ) THEN
                        CALL SETMSG ( 'Integer error (#) in line #.' )
                        CALL ERRCH  ( '#', ERROR                     )
                        CALL ERRINT ( '#', NREC                      )
                        CALL SIGERR ( 'SPICE(BADINTEGER)'            )
                        CALL CHKOUT ( 'RDFFPL'                       )
                        RETURN
                     END IF

                  END DO

                  PLATES(1,NP) = ARI(1)
                  PLATES(2,NP) = ARI(2)
                  PLATES(3,NP) = ARI(3)

               ELSE

C
C                 This block executes if the first characters of 
C                 a type 3 data lacks a 'v' or 'f'.
C
                  CALL SETMSG ( 'Bad data line at record #.' )
                  CALL ERRINT ( '#', NREC                    )
                  CALL SIGERR ( 'SPICE(BADDATALINE)'         )
                  CALL CHKOUT ( 'RDFFPL'                     )
                  RETURN

               END IF

               NREC = NREC + 1

            END IF

         END DO

C
C        All type 3 data read. We should have a non zero NP and NV,
C        if not, something failed.
C
         IF ( (NP .EQ. 0) .OR. (NV .EQ. 0) ) THEN
            CALL SETMSG ( 'Read error, type 3 data file.'
     .                 // ' Num plates = #1, num verts = #2.'
     .                 // 'Both should be non-zero'         )
            CALL ERRINT ( '#1', NP                          )
            CALL ERRINT ( '#2', NV                          )
            CALL SIGERR ( 'SPICE(DATAREADFAILED)'           )
            CALL CHKOUT ( 'RDFFPL'                          )
            RETURN
         END IF





      ELSE IF ( PLTTYP .EQ. 4 ) THEN
C
C        We have a Rosetta Osiris style ".ver" file.
C
C        Get the vertex and plate count from the first non-blank
C        input line.
C
         CALL RDNBL ( INFILE, ILINE, EOF )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'RDFFPL' )
            RETURN
         END IF

         CALL LPARSM ( ILINE, ' ,', 2, NTOK, TOKENS )

         IF ( NTOK .NE. 2 ) THEN
         
            CALL SETMSG ( 'Vertex and plate count were expected '
     .      //            'on first line of input file. Line was '
     .      //            '#.'                                    )
            CALL ERRCH  ( '#', ILINE                              )
            CALL SIGERR ( 'SPICE(UNRECOGNIZEDFORMAT)'             )
            CALL CHKOUT ( 'RDFFPL'                                )
            RETURN

         END IF

         CALL PRSINT ( TOKENS(1), NV )
         CALL PRSINT ( TOKENS(2), NP )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'RDFFPL' )
            RETURN
         END IF


         IF ( NV .GT. MAXVRT ) THEN

            CALL SETMSG ( 'Number of vertices # exceeds limit #.' )
            CALL ERRINT ( '#',  NV                                )
            CALL ERRINT ( '#',  MAXVRT                            )
            CALL SIGERR ( 'SPICE(TOOMANYVERTICES)'                )
            CALL CHKOUT ( 'RDFFPL'                                )
            RETURN

         END IF

         IF ( NP .GT. MAXPLT ) THEN

            CALL SETMSG ( 'Number of plates # exceeds limit #.' )
            CALL ERRINT ( '#', NP                               )
            CALL ERRINT ( '#', MAXPLT                           )
            CALL SIGERR ( 'SPICE(TOOMANYPLATES)'                )
            CALL CHKOUT ( 'RDFFPL'                              )
            RETURN

         END IF
 
C
C        Read the vertex data and store it in the output vertex array.
C
         DO I = 1, NV
 
            CALL RDNBL ( INFILE, ILINE, EOF )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'RDFFPL' )
               RETURN
            END IF
 
            IF ( EOF ) THEN

               CALL SETMSG ( 'Expected to find # vertices in input '
     .         //            'file # but ran out of data after reading '
     .         //            '# lines of vertex data.'                 )
               CALL ERRINT ( '#', NV                                   )
               CALL ERRCH  ( '#', INFILE                               )
               CALL ERRINT ( '#', I                                    )
               CALL SIGERR ( 'SPICE(FILETRUNCATED)'                    )
               CALL CHKOUT ( 'RDFFPL'                                  )
               RETURN

            END IF
 
C
C           Parse the input line; we expect it to contain the 
C           components of the Ith vertex.
C
            CALL LPARSM ( ILINE, ' ,', 3, NTOK, TOKENS )

            IF ( NTOK .NE. 3 ) THEN
         
               CALL SETMSG ( 'Three vertex components were expected '
     .         //            'on current line of input file. Line was '
     .         //            '#.'                                    )
               CALL ERRCH  ( '#', ILINE                              )
               CALL SIGERR ( 'SPICE(UNRECOGNIZEDFORMAT)'             )
               CALL CHKOUT ( 'RDFFPL'                                )
               RETURN

            END IF

         
            CALL PRSDP ( TOKENS(1), VRTCES(1,I) )
            CALL PRSDP ( TOKENS(2), VRTCES(2,I) )
            CALL PRSDP ( TOKENS(3), VRTCES(3,I) )

            IF ( FAILED() ) THEN
               CALL CHKOUT( 'RDFFPL' )
               RETURN
            END IF

         END DO
 
C
C        Read the plate indices. Discard the first
C        and every second non-blank line of plate data; these
C        lines all contain the same vertex count (3).
C
         J = 0

         DO I = 1,  2 * NP

            CALL RDNBL ( INFILE, ILINE, EOF )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'RDFFPL' )
               RETURN
            END IF

            IF ( EOF ) THEN

               CALL SETMSG ( 'Expected to find # plates in input '
     .         //            'file # but ran out of data after reading '
     .         //            '# lines of plate data.'                  )
               CALL ERRINT ( '#', NP                                   )
               CALL ERRCH  ( '#', INFILE                               )
               CALL ERRINT ( '#', I                                    )
               CALL SIGERR ( 'SPICE(FILETRUNCATED)'                    )
               CALL CHKOUT ( 'RDFFPL'                                  )
               RETURN

            END IF

            IF ( EVEN(I) ) THEN
C
C              This line should contain the components of the Jth plate.
C
               J = J + 1

C
C              Parse the current line and store the plate components.
C
               CALL LPARSM ( ILINE, ' ,', 3, NTOK, TOKENS )

               IF ( NTOK .NE. 3 ) THEN
         
                  CALL SETMSG ( 'Three plate components were expected '
     .            //            'on current line of input file. Line '
     .            //            'was #.'                              )
                  CALL ERRCH  ( '#', ILINE                            )
                  CALL SIGERR ( 'SPICE(UNRECOGNIZEDFORMAT)'           )
                  CALL CHKOUT ( 'RDFFPL'                              )
                  RETURN

               END IF
         
               CALL PRSINT ( TOKENS(1), PLATES(1,J) )
               CALL PRSINT ( TOKENS(2), PLATES(2,J) )
               CALL PRSINT ( TOKENS(3), PLATES(3,J) )

               IF ( FAILED() ) THEN
                  CALL CHKOUT( 'RDFFPL' )
                  RETURN
               END IF

            END IF

         END DO


      ELSE

         CALL SETMSG ( 'Unkown plate data type: #.' )
         CALL ERRINT( '#', PLTTYP                   )
         CALL SIGERR ( 'SPICE(BADDATATYPE)'         )
         CALL CHKOUT ( 'RDFFPL'                     )
         RETURN

      END IF

      CALL TOSTDO ( '...Done reading plate model input file.' )
      CALL TOSTDO ( ' ' )
 
      CALL CHKOUT ( 'RDFFPL' )
 
      RETURN
      END




