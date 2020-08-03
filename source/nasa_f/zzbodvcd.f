C$Procedure ZZBODVCD ( Get d.p. kernel variable for body, with bypass )
 
      SUBROUTINE ZZBODVCD ( BODYID, ITEM, MAXN, VARCTR, N, VALUES )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Fetch from the kernel pool the double precision values of an item
C     associated with a body, where the body is specified by an integer
C     ID code. Perform lookup only if kernel variable or kernel pool
C     state has changed; otherwise, bypass kernel pool lookup and
C     return input values.
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
C     KERNEL
C     NAIF_IDS
C
C$ Keywords
C
C     CONSTANTS
C
C$ Declarations
 
      IMPLICIT NONE
      INCLUDE 'zzctr.inc'

      INTEGER               BODYID
      CHARACTER*(*)         ITEM
      INTEGER               MAXN
      INTEGER               VARCTR ( CTRSIZ )
      INTEGER               N
      DOUBLE PRECISION      VALUES ( * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   Body ID code.
C     ITEM       I   Item for which values are desired.
C     MAXN       I   Maximum number of values that may be returned.
C     VARCTR    I-O  POOL state counter saved by the caller.
C     N         I-O  Number of values returned.
C     VALUES    I-O  Array of double precision values.
C
C$ Detailed_Input
C
C     BODYID     is the NAIF integer ID code for a body of interest.
C                For example, if the body is the earth, the code is
C                399.
C
C
C     ITEM       is the item to be returned. Together, the NAIF ID
C                code of the body and the item name combine to form a
C                kernel variable name, e.g.,
C
C                      'BODY599_RADII'     
C                      'BODY401_POLE_RA' 
C
C                The values associated with the kernel variable having
C                the name constructed as shown are sought. Below
C                we'll take the shortcut of calling this kernel variable
C                the "requested kernel variable."
C
C                Note that ITEM *is* case-sensitive. This attribute
C                is inherited from the case-sensitivity of kernel
C                variable names.
C
C
C     MAXN       is the maximum number of values that may be returned.
C                The output array VALUES must be declared with size at
C                least MAXN.  It's an error to supply an output array
C                that is too small to hold all of the values associated
C                with the requested kernel variable.
C
C
C     VARCTR     is, on input, a counter used to represent the kernel
C                pool state at the time the kernel variable designated
C                by BODYID and ITEM was set.
C
C                Note that a distinct actual argument corresponding to
C                VARCTR is required for each kernel variable to be
C                updated by this routine.
C
C
C     N          is, on input, the dimension of the saved kernel
C                variable designated by BODYID and ITEM.
C
C
C     VALUES     is, on input, the set of values associated with the
C                kernel variable designated by BODYID and ITEM.
C
C$ Detailed_Output
C
C     VARCTR     is, on output, a counter used to represent the current
C                kernel pool state.
C
C     N          is, on output, the dimension of the kernel variable
C                designated by BODYID and ITEM.
C
C     VALUES     is, on output, the set of values associated with the
C                kernel variable designated by BODYID and ITEM.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the requested kernel variable is not found in the kernel
C         pool, the error SPICE(KERNELVARNOTFOUND) is signaled.
C
C     2)  If the requested kernel variable is found but the associated
C         values aren't numeric, an error is signaled by a routine in
C         the call tree of this routine.
C
C     3)  The output array VALUES must be declared with sufficient size
C         to contain all of the values associated with the requested
C         kernel variable.  If the dimension of VALUES indicated by
C         MAXN is too small to contain the requested values, an error
C         is signaled by a routine in the call tree of this routine.
C
C     4)  If the input dimension MAXN indicates there is more room in
C         VALUES than there really is---for example, if MAXN is 10 but
C         values is declared with dimension 5---and the dimension of
C         the requested kernel variable is larger than the actual
C         dimension of VALUES, then this routine may overwrite memory.
C         The results are unpredictable.
C
C     5)  If a lookup of the requested variable fails, the argument
C         N is set to zero.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is meant to enable efficient use of BODVCD. This
C     routine calls BODVCD only if the requested kernel variable values
C     may have changed from their previous values.
C
C     This routine simplifies looking up PCK kernel variables by
C     constructing names of requested kernel variables and by
C     performing error checking.
C
C     This routine is intended for use in cases where the maximum
C     number of values that may be returned is known at compile time.
C     The caller fetches all of the values associated with the
C     specified kernel variable via a single call to this routine. If
C     the number of values to be fetched cannot be known until run
C     time, the lower-level routine GDPOOL (an entry point of POOL)
C     should be used instead. GDPOOL supports fetching arbitrary
C     amounts of data in multiple "chunks."
C
C     This routine is intended for use in cases where the requested
C     kernel variable is expected to be present in the kernel pool. If
C     the variable is not found or has the wrong data type, this
C     routine signals an error. In cases where it is appropriate to
C     indicate absence of an expected kernel variable by returning a
C     boolean "found flag" with the value .FALSE., again the routine
C     GDPOOL should be used.
C
C$ Examples
C
C     1)  The call below can be used to fetch initial values of the
C         radii for the earth (body 399).
C
C         When the kernel variable 
C
C            BODY399_RADII
C
C         is present in the kernel pool---normally because a PCK
C         defining this variable has been loaded---the code fragment    
C      
C            CALL ZZCTRUIN ( CTR1 )
C
C            BODYID = 399
C            ITEM   = 'RADII'
C
C            CALL ZZBODVCD ( BODYID, ITEM, MAXN, CTR1, N, RADII )
C
C         returns the dimension and values associated with the variable
C         'BODY399_RADII', for example,
C
C            DIM      = 3
C            RADII(1) = 6378.140
C            RADII(2) = 6378.140
C            RADII(3) = 6356.755
C
C         The call to ZZCTRUIN initializes CTR1, which forces ZZBODVCD
C         to update the output values
C
C
C     2)  After the call in example (1) has been made, the variable
C
C            BODY399_RADII
C
C         can be set to its current value by the call
C
C            CALL ZZBODVCD ( BODYID, ITEM, MAXN, CTR1, N, RADII )
C
C         If the kernel pool has not been modified since the previous
C         call, the arguments
C
C            CTR1
C            N
C            RADII
C
C         are unchanged on output.
C
C
C     3)  After the call in example (2) has been made, the variable
C
C            BODY499_POLE_RA
C
C         can be set to its current value by the code fragment
C
C            CALL ZZCTRUIN ( CTR2 )
C
C            BODY   = 499
C            ITEM   = 'POLE_RA'
C
C            CALL ZZBODVCD ( BODYID, ITEM, MAXN, CTR2, N, RA )
C
C         On output, the arguments
C
C            N
C            RA
C
C         will be updated to (based on PCK pck00010.tpc):
C
C            3
C            499
C            'POLE_RA'
C            317.68143   -0.1061      0.  
C
C         CTR2 will be updated if the kernel pool contents have
C         changed.
C        
C
C     4) The call 
C
C            BODY   = 399
C            ITEM   = 'radii'
C
C            CALL ZZBODVCD ( BODYID, ITEM, MAXN, CTR1, N, VALUES )
C
C        usually will cause a SPICE(KERNELVARNOTFOUND) error to be
C        signaled, because this call will attempt to look up the
C        values associated with a kernel variable of the name
C
C           'BODY399_radii'
C
C        Since kernel variable names are case sensitive, this
C        name is not considered to match the name
C
C           'BODY399_RADII'
C
C        which normally would be present after a text PCK
C        containing data for all planets and satellites has 
C        been loaded.
C 
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C     W.L. Taber      (JPL)
C     I.M. Underwood  (JPL)
C
C$ Version
C
C-    SPICELIB Version 2.0.0, 20-OCT-2015 (NJB) (BVS) (WLT) (IMU)
C
C-&
 
C$ Index_Entries
C
C     fetch constants for a body from the kernel pool
C     physical constants for a body
C
C-&
 
C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C
      LOGICAL               UPDATE
 
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZBODVCD' )
            
C
C     See whether the kernel pool state has changed since the 
C     user counter was set. Update the user counter if so.
C
      CALL ZZPCTRCK ( VARCTR, UPDATE )

C
C     If the pool was updated, or if we're looking at a new variable,
C     update the kernel variable values, size, and found flag. 
C     Otherwise do nothing.
C
      IF ( UPDATE ) THEN

         CALL BODVCD ( BODYID, ITEM, MAXN, N, VALUES )

         IF ( FAILED() ) THEN

            N = 0

         END IF

      END IF

      CALL CHKOUT ( 'ZZBODVCD' )
      RETURN
      END
