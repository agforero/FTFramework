C$Procedure ZZSRFKER ( Surface translation, process kernel update )
 
      SUBROUTINE ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                      EXTKER, NKVAR,  SNMHLS, SNMPOL, 
     .                      SNMIDX, SIDHLS, SIDPOL, SIDIDX )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Update ZZSRFTRN's name-based and ID-based data structure arrays
C     using the contents of kernel variables that define the surface
C     name/ID mapping.
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
C     None.
C
C$ Keywords
C
C     UTILITY
C
C$ Declarations
 
      IMPLICIT NONE
 
      INCLUDE              'srftrn.inc'
      INCLUDE              'zzsrftrn.inc'
 
      CHARACTER*(SFNMLN)    KERNAM (          * )
      CHARACTER*(SFNMLN)    NORNAM (          * )
      INTEGER               KERSID (          * )      
      INTEGER               KERBID (          * )
      LOGICAL               EXTKER
      INTEGER               NKVAR
      INTEGER               SNMHLS (          * )
      INTEGER               SNMPOL ( LBSNGL : * )
      INTEGER               SNMIDX (          * )
      INTEGER               SIDHLS (          * )
      INTEGER               SIDPOL ( LBSNGL : * )
      INTEGER               SIDIDX (          * )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     KERNAM     O   Array of surface names from kernel pool.
C     NORNAM     O   Array of normalized surface names.
C     KERSID     O   Array of surface ID codes from kernel pool.
C     KERBID     O   Array of body ID codes from kernel pool.
C     EXTKER     O   Logical flag indicating kernel data are present.
C     NKVAR      O   Number of surface name, code, body tuples.
C     SNMHLS     O   Surface name-based hash head node pointer list
C     SNMPOL     O   Surface name-based hash node collision list
C     SNMIDX     O   Surface name-based hash index storage array 
C     SIDHLS     O   Surface ID-based hash head node pointer list
C     SIDPOL     O   Surface ID-based hash node collision list
C     SIDIDX     O   Surface ID-based hash index storage array
C     LBSNGL     P   Lower bound of hash pool arrays
C     SFNMLN     P   Maximum length of surface name strings
C
C$ Detailed_Input
C
C     None.
C
C$ Detailed_Output
C
C     KERNAM    is an array containing surface names from kernel pool
C               assignments to the kernel variable NAIF_SURFACE_NAME.
C
C               Array elements from masked assignments are not included
C               in KERNAM.
C
C     NORNAM    is an array parallel to KERNAM containing normalized
C               names. The Ith element of NORNAM is obtained from the
C               Ith element of KERNAM by conversion to uppercase, 
C               left-justification, and compression of consecutive
C               embedded blanks to a single blank.
C
C     KERSID    is an array containing surface names from kernel pool
C               assignments to the kernel variable NAIF_SURFACE_CODE.
C               The Ith element of KERSID is the code associated with
C               the Ith element of KERNAM.
C
C     KERBID    is an array containing surface names from kernel pool
C               assignments to the kernel variable NAIF_SURFACE_BODY.
C               The Ith element of KERBID is the code associated with
C               the Ith element of KERNAM.
C
C     EXTKER    is a logical flag indicating whether kernel data
C               defining a surface name/ID mapping are available.
C               EXTKER is set to .TRUE. if the data are present and is
C               .FALSE. otherwise.
C
C     NKVAR     is the count of names in the array KERNAM; the
C               arrays NORNAM, KERSID, and KERBID also contain NKVAR
C               entries.
C
C     SNMHLS
C     SNMPOL    are the surface name-based hash head node pointer and
C               collision lists. Together with the arrays SNMIDX,
C               NORNAM, KERBID, and KERSID, they enable mapping pairs
C               of normalized surface names and body ID codes to
C               surface ID codes.
C
C     SNMIDX    is the surface name-based hash index storage array. It
C               maps nodes in the name collision list to entries in the 
C               parallel NORNAM, KERBID, and KERSID arrays.
C
C     SIDHLS
C     SIDPOL    are the surface ID-based hash head node pointer and  
C               collision lists. Together with the arrays SIDIDX,
C               KERSID, KERBID, and KERNAM, they enable mapping pairs of
C               surface ID codes and body ID codes to surface names.
C
C     SIDIDX    is the surface ID-based hash index storage array. It
C               maps nodes in the ID collision list to entries in the 
C               parallel KERSID, KERBID, and KERNAM arrays.
C 
C$ Parameters
C
C     LBSNGL    is the lower bound of the hashes' collision list array.
C
C     SFNMLN    is the maximum length of a surface name. Defined in the
C               include file 'srftrn.inc'.
C
C$ Exceptions
C
C     1)  If an error occurs while fetching kernel variables, the error
C         will be signaled by a routine in the call tree of this
C         routine.
C
C     2)  All three of the kernel variables defining the surface name/ID
C         map must be present in the kernel pool, or all must be absent.
C         If this condition is not met, the error SPICE(BADSURFACEMAP)
C         will be signaled.
C
C     3)  All three of the kernel variables defining the surface name/ID
C         map must have the correct data type; if not, the error
C         SPICE(BADVARIABLETYPE) will be signaled.
C
C     4)  If any of the kernel variables defining the surface name/ID
C         map have size exceeding MXNSRF, the error
C         SPICE(TOOMANYSURFACES) will be signaled.
C
C     5)  All three of the kernel variables defining the surface name/ID
C         map must have the same number of elements; if not, the error
C         SPICE(ARRAYSIZEMISMATCH) will be signaled.
C
C     6)  If any surface name is blank, the error
C         SPICE(BLANKNAMEASSIGNED) will be signaled.
C
C     7)  Any error occurring during initialization of the hash data
C         structures will be signaled by a routine in the call tree
C         of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine manages updates of the surface name/ID mapping
C     subsystem data structures. The outputs of this routine contain
C     mapping information provided by kernel variables, as well as
C     associated bookkeeping information.
C
C     On the first pass through this routine, the routine sets a watch
C     on the kernel variables that define the surface name/ID mapping.
C     This routine must be called after any change to these variables.
C     
C$ Examples
C
C     See the routine ZZSRFTRN.
C
C$ Restrictions
C
C     1)  This routine is intended only for use by ZZSRFTRN and
C         ZZSRFKER.
C
C     2)  All output hash arrays must be declared with the same
C         dimension which is greater than or equal to MAXVAL.
C
C     3)  The order of mappings in the input arrays determines the
C         priority, with the mapping with the lowest priority being the
C         first and the mapping with the highest priority being the
C         last.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     B.V. Semenov       (JPL)
C     M.J. Spencer       (JPL)
C     W.L. Taber         (JPL)
C     F.S. Turner        (JPL)
C     E.D. Wright        (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 04-DEC-2015 (NJB) (BVS) (EDW)
C
C-&
 


C
C     SPICELIB functions
C      
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      CHARACTER*(*)         KVSFBD
      PARAMETER           ( KVSFBD = 'NAIF_SURFACE_BODY' )

      CHARACTER*(*)         KVSFCD
      PARAMETER           ( KVSFCD = 'NAIF_SURFACE_CODE' )

      CHARACTER*(*)         KVSFNM
      PARAMETER           ( KVSFNM = 'NAIF_SURFACE_NAME' )

      INTEGER               KVNMLN
      PARAMETER           ( KVNMLN = 32 )

      INTEGER               NNAMES
      PARAMETER           ( NNAMES = 3 )

C
C     Local variables
C
      CHARACTER*(1)         BDTYPE
      CHARACTER*(1)         CDTYPE
      CHARACTER*(1)         NDTYPE
      CHARACTER*(KVNMLN)    NAMES  ( NNAMES )

      INTEGER               I
      INTEGER               NBODY
      INTEGER               NCODE
      INTEGER               NNAME

      LOGICAL               PASS1
      LOGICAL               FNDBOD
      LOGICAL               FNDCDE
      LOGICAL               FNDNAM


C
C     Saved variables
C
      SAVE                  PASS1

C
C     Initial values
C     
      DATA                  PASS1 / .TRUE. /
      DATA                  NAMES / KVSFBD, KVSFCD, KVSFNM /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSRFKER' )

C
C     The primary functions performed inline in this routine are
C
C        - Setting a watch on the mapping kernel variables
C
C        - Fetching the mapping kernel variables' values
C
C        - Performing error checks on the kernel variables that define
C          the surface name/ID mapping. Initialization of data
C          structures is delegated to ZZSRFINI.
C
C
      IF ( PASS1 ) THEN
C
C        Set watch on kernel variables used for the surface mapping.
C
         CALL SWPOOL ( 'ZZSRFTRN', NNAMES, NAMES )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSRFKER' )
            RETURN
         END IF

         PASS1 = .FALSE.

      END IF

C
C     Indicate that no data are available until we find out
C     otherwise.
C
      EXTKER = .FALSE.
      NKVAR  =  0

C
C     Fetch attributes of the surface mapping kernel variables.
C
      CALL DTPOOL ( KVSFNM, FNDNAM, NNAME, NDTYPE )
      CALL DTPOOL ( KVSFCD, FNDCDE, NCODE, CDTYPE )
      CALL DTPOOL ( KVSFBD, FNDBOD, NBODY, BDTYPE )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSRFKER' )
         RETURN
      END IF

C
C     The variables must all be present or all be absent.
C
      IF ( ( FNDCDE .NEQV. FNDNAM ) .OR. ( FNDBOD .NEQV. FNDNAM ) ) THEN

         CALL SETMSG ( 'Surface mapping kernel variables are in an '
     .   //            'inconsistent state. # was #; # was #; # '
     .   //            'was #.'                                      )

         CALL ERRCH  ( '#', KVSFNM )

         IF ( FNDNAM ) THEN
            CALL ERRCH ( '#', 'found' )
         ELSE 
            CALL ERRCH ( '#', 'not found' )
         END IF


         CALL ERRCH  ( '#', KVSFCD )

         IF ( FNDCDE ) THEN
            CALL ERRCH ( '#', 'found' )
         ELSE 
            CALL ERRCH ( '#', 'not found' )
         END IF


         CALL ERRCH  ( '#', KVSFBD )

         IF ( FNDBOD ) THEN
            CALL ERRCH ( '#', 'found' )
         ELSE 
            CALL ERRCH ( '#', 'not found' )
         END IF

         CALL SIGERR ( 'SPICE(BADSURFACEMAP)' )
         CALL CHKOUT ( 'ZZSRFKER'             )
         RETURN

      END IF

C
C     If the variables are not present, leave now.
C
      EXTKER = FNDNAM .AND. FNDCDE .AND. FNDBOD

      IF ( .NOT. EXTKER ) THEN
         CALL CHKOUT ( 'ZZSRFKER' )
         RETURN
      END IF

C
C     Make sure the kernel variables aren't larger than our arrays.
C     Also make sure the variable have matching dimensions.
C
C     Check variable types.
C
      IF (     ( NDTYPE .NE. 'C' )
     .    .OR. ( CDTYPE .NE. 'N' )
     .    .OR. ( BDTYPE .NE. 'N' ) ) THEN

         CALL SETMSG ( 'Surface mapping kernel variable types '
     .   //            'are: # = #; # = #; # = #. These types '
     .   //            'must be, respectively, ''C'', ''N'', '
     .   //            '''N''.'                                )
         CALL ERRCH  ( '#', KVSFNM                             )
         CALL ERRCH  ( '#', NDTYPE                             )
         CALL ERRCH  ( '#', KVSFCD                             )
         CALL ERRCH  ( '#', CDTYPE                             )
         CALL ERRCH  ( '#', KVSFBD                             )
         CALL ERRCH  ( '#', BDTYPE                             )
         CALL SIGERR ( 'SPICE(BADVARIABLETYPE)'                )
         CALL CHKOUT ( 'ZZSRFKER'                              )
         RETURN

      END IF

C
C     Check variable dimensions.
C     
      IF (     ( NNAME .GT. MXNSRF )
     .    .OR. ( NCODE .GT. MXNSRF )
     .    .OR. ( NBODY .GT. MXNSRF ) ) THEN

         CALL SETMSG ( 'Surface mapping kernel variable sizes '
     .   //            'are: # = #; # = #; # = #. Maximum '
     .   //            'allowed size is #.'                    )
         CALL ERRCH  ( '#', KVSFNM                             )
         CALL ERRINT ( '#', NNAME                              )
         CALL ERRCH  ( '#', KVSFCD                             )
         CALL ERRINT ( '#', NCODE                              )
         CALL ERRCH  ( '#', KVSFBD                             )
         CALL ERRINT ( '#', NBODY                              )
         CALL ERRINT ( '#', MXNSRF                             )
         CALL SIGERR ( 'SPICE(TOOMANYSURFACES)'                )
         CALL CHKOUT ( 'ZZSRFKER'                              )
         RETURN

      END IF

      IF (  ( NCODE .NE. NNAME ) .OR. ( NBODY .NE. NNAME )  ) THEN

         CALL SETMSG ( 'Surface variable sizes do not match. Size '
     .   //            'of # is #; size of # is #; size of # is #.'  ) 
         CALL ERRCH  ( '#', KVSFNM                                   )
         CALL ERRINT ( '#', NNAME                                    )
         CALL ERRCH  ( '#', KVSFCD                                   )
         CALL ERRINT ( '#', NCODE                                    )
         CALL ERRCH  ( '#', KVSFBD                                   )
         CALL ERRINT ( '#', NBODY                                    )
         CALL SIGERR ( 'SPICE(ARRAYSIZEMISMATCH)'                    )
         CALL CHKOUT ( 'ZZSRFKER'                                    )
         RETURN

      END IF

C
C     Fetch mapping variables.
C
C     Note that we'll check the variable sizes below.
C
      CALL GCPOOL ( KVSFNM, 1, MXNSRF, NNAME, KERNAM, FNDNAM )
      CALL GIPOOL ( KVSFCD, 1, MXNSRF, NCODE, KERSID, FNDCDE )
      CALL GIPOOL ( KVSFBD, 1, MXNSRF, NBODY, KERBID, FNDBOD )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSRFKER' )
         RETURN
      END IF


      NKVAR = NNAME

C
C     Produce normalized name array. Check for blank names 
C     as we go.
C
      DO I = 1, NKVAR

         IF ( KERNAM(I) .EQ. ' ' ) THEN

            CALL SETMSG ( 'An attempt to assign the code, #, to '
     .      //            'a blank string was made.  Check loaded '
     .      //            'text kernels for a blank string in '
     .      //            'the NAIF_SURFACE_NAME array.'            )
            CALL ERRINT ( '#', I                                    )
            CALL SIGERR ( 'SPICE(BLANKNAMEASSIGNED)'                )
            CALL CHKOUT ( 'ZZSRFKER'                                )
            RETURN

         END IF

         CALL LJUCRS ( 1, KERNAM(I), NORNAM(I) )

      END DO

C
C     Initialize hash data structures.
C     
      CALL ZZSRFINI ( NORNAM, KERSID, KERBID, NKVAR, 
     .                NROOM,  SNMHLS, SNMPOL, SNMIDX,
     .                SIDHLS, SIDPOL, SIDIDX         )
 
      CALL CHKOUT ( 'ZZSRFKER' )
      RETURN
      END
