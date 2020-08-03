C$Procedure      BODMAT ( Return transformation matrix for a body )
 
      SUBROUTINE BODMAT ( BODY, ET, TIPM )
 
C$ Abstract
C
C     Return the J2000 to body Equator and Prime Meridian coordinate
C     transformation matrix for a specified body.
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
C     PCK
C     NAIF_IDS
C     TIME
C
C$ Keywords
C
C     CONSTANTS
C
C$ Declarations

      IMPLICIT NONE
 
      INTEGER               BODY
      DOUBLE PRECISION      ET
      DOUBLE PRECISION      TIPM   ( 3,3 )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODY       I   ID code of body.
C     ET         I   Epoch of transformation.
C     TIPM       O   Transformation from Inertial to PM for BODY at ET.
C
C$ Detailed_Input
C
C     BODY        is the integer ID code of the body for which the
C                 transformation is requested. Bodies are numbered
C                 according to the standard NAIF numbering scheme.
C
C     ET          is the epoch at which the transformation is
C                 requested. (This is typically the epoch of
C                 observation minus the one-way light time from
C                 the observer to the body at the epoch of
C                 observation.)
C
C$ Detailed_Output
C
C     TIPM        is the transformation matrix from Inertial to body
C                 Equator and Prime Meridian.  The X axis of the PM
C                 system is directed to the intersection of the
C                 equator and prime meridian. The Z axis points north.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If data required to define the body-fixed frame associated
C        with BODY are not found in the binary PCK system or the kernel
C        pool, the error SPICE(FRAMEDATANOTFOUND) is signaled. In
C        the case of IAU style body-fixed frames, the absence of
C        prime meridian polynomial data (which are required) is used
C        as an indicator of missing data.
C
C     2) If the test for exception (1) passes, but in fact requested
C        data are not available in the kernel pool, the error will be
C        signaled by routines in the call tree of this routine.
C
C     3) If the kernel pool does not contain all of the data required
C        to define the number of nutation precession angles
C        corresponding to the available nutation precession
C        coefficients, the error SPICE(INSUFFICIENTANGLES) is
C        signaled.
C
C     4) If the reference frame REF is not recognized, a routine
C        called by BODMAT will diagnose the condition and invoke the
C        SPICE error handling system.
C
C     5) If the specified body code BODY is not recognized, the
C        error is diagnosed by a routine called by BODMAT.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is related to the more general routine TIPBOD
C     which returns a matrix that transforms vectors from a
C     specified inertial reference frame to body equator and
C     prime meridian coordinates.  TIPBOD accepts an input argument
C     REF that allows the caller to specify an inertial reference
C     frame.
C
C     The transformation represented by BODMAT's output argument TIPM
C     is defined as follows:
C
C        TIPM = [W] [DELTA] [PHI]
C                 3        1     3
C
C     If there exists high-precision binary PCK kernel information
C     for the body at the requested time, these angles, W, DELTA
C     and PHI are computed directly from that file.  The most
C     recently loaded binary PCK file has first priority followed
C     by previously loaded binary PCK files in backward time order.
C     If no binary PCK file has been loaded, the text P_constants
C     kernel file is used.
C
C     If there is only text PCK kernel information, it is
C     expressed in terms of RA, DEC and W (same W as above), where
C
C        RA    = PHI - HALFPI()
C        DEC   = HALFPI() - DELTA
C
C     RA, DEC, and W are defined as follows in the text PCK file:
C
C           RA  = RA0  + RA1*T  + RA2*T*T   + a  sin theta
C                                              i          i
C
C           DEC = DEC0 + DEC1*T + DEC2*T*T  + d  cos theta
C                                              i          i
C
C           W   = W0   + W1*d   + W2*d*d    + w  sin theta
C                                              i          i
C
C     where:
C
C           d = days past J2000.
C
C           T = Julian centuries past J2000.
C
C           a , d , and w  arrays apply to satellites only.
C            i   i       i
C
C           theta  = THETA0 * THETA1*T are specific to each planet.
C                i
C
C     These angles -- typically nodal rates -- vary in number and
C     definition from one planetary system to the next.
C
C$ Examples
C
C     In the following code fragment, BODMAT is used to rotate
C     the position vector (POS) from a target body (BODY) to a
C     spacecraft from inertial coordinates to body-fixed coordinates
C     at a specific epoch (ET), in order to compute the planetocentric
C     longitude (PCLONG) of the spacecraft.
C
C        CALL BODMAT ( BODY, ET, TIPM )
C        CALL MXV    ( TIPM, POS, POS )
C        CALL RECLAT ( POS, RADIUS, PCLONG, LAT )
C
C     To compute the equivalent planetographic longitude (PGLONG),
C     it is necessary to know the direction of rotation of the target
C     body, as shown below.
C
C        CALL BODVCD ( BODY, 'PM', 3, DIM, VALUES )
C
C        IF ( VALUES(2) .GT. 0.D0 ) THEN
C           PGLONG = PCLONG
C        ELSE
C           PGLONG = TWOPI() - PCLONG
C        END IF
C
C     Note that the items necessary to compute the transformation
C     TIPM must have been loaded into the kernel pool (by one or more
C     previous calls to FURNSH).
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     1)  Refer to the NAIF_IDS required reading file for a complete
C         list of the NAIF integer ID codes for bodies.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C     W.L. Taber      (JPL)
C     I.M. Underwood  (JPL)
C     K.S. Zukor      (JPL)
C
C$ Version
C
C-    SPICELIB Version 4.2.0, 27-JUL-2016 (BVS)
C
C        Updated to use the 3x3 top-left corner of the 6x6 matrix
C        returned by TISBOD instead of fetching kernel data and doing
C        computations in-line.
C
C-    SPICELIB Version 4.1.1, 01-FEB-2008 (NJB)
C
C        The routine was updated to improve the error messages created
C        when required PCK data are not found. Now in most cases the
C        messages are created locally rather than by the kernel pool
C        access routines. In particular missing binary PCK data will
C        be indicated with a reasonable error message.
C
C-    SPICELIB Version 4.1.0, 25-AUG-2005 (NJB)
C
C        Updated to remove non-standard use of duplicate arguments
C        in MXM call.
C
C         Calls to ZZBODVCD have been replaced with calls to 
C         BODVCD.
C
C-     SPICELIB Version 4.0.0, 12-FEB-2004 (NJB)
C
C         Code has been updated to support satellite ID codes in the
C         range 10000 to 99999 and to allow nutation precession angles
C         to be associated with any object. 
C
C         Implementation changes were made to improve robustness
C         of the code.
C
C-     SPICELIB Version 3.2.0, 22-MAR-1995 (KSZ)
C
C        Gets TSIPM matrix from PCKMAT (instead of Euler angles
C        from PCKEUL.)
C
C-     SPICELIB Version 3.0.0, 10-MAR-1994 (KSZ)
C
C        Ability to get Euler angles from binary PCK file added.
C        This uses the new routine PCKEUL.
C
C-     SPICELIB Version 2.0.1, 10-MAR-1992 (WLT)
C
C         Comment section for permuted index source lines was added
C         following the header.
C
C-     SPICELIB Version 2.0.0, 04-SEP-1991 (NJB)
C
C         Updated to handle P_constants referenced to different epochs
C         and inertial reference frames.
C
C         The header was updated to specify that the inertial reference
C         frame used by BODMAT is restricted to be J2000.
C
C-    SPICELIB Version 1.0.0, 31-JAN-1990 (WLT) (IMU)
C
C-&
 
C$ Index_Entries
C
C     fetch transformation matrix for a body
C     transformation from j2000 position to bodyfixed
C     transformation from j2000 to bodyfixed coordinates
C
C-&
 
 
C$ Revisions
C
C-    SPICELIB Version 4.2.0, 02-MAR-2016 (BVS)
C
C        Updated to use the 3x3 top-left corner of the 6x6 matrix
C        returned by TISBOD instead of fetching kernel data and doing
C        computations in-line.
C
C-    SPICELIB Version 4.1.0, 25-AUG-2005 (NJB)
C
C        Updated to remove non-standard use of duplicate arguments
C        in MXM call.
C
C         Calls to ZZBODVCD have been replaced with calls to 
C         BODVCD.
C
C-     SPICELIB Version 4.0.0, 12-FEB-2004 (NJB)
C
C         Code has been updated to support satellite ID codes in the
C         range 10000 to 99999 and to allow nutation precession angles
C         to be associated with any object.
C
C         Calls to deprecated kernel pool access routine RTPOOL 
C         were replaced by calls to GDPOOL.
C
C         Calls to BODVAR have been replaced with calls to 
C         ZZBODVCD.
C
C-     SPICELIB Version 3.2.0, 22-MAR-1995 (KSZ)
C
C        BODMAT now get the TSIPM matrix from PCKMAT, and
C        unpacks TIPM from it.  Also the calculated but unused
C        variable LAMBDA was removed.
C
C-     SPICELIB Version 3.0.0, 10-MAR-1994 (KSZ)
C
C        BODMAT now uses new software to check for the
C        existence of binary PCK files, search the for
C        data corresponding to the requested body and time,
C        and return the appropriate Euler angles, using the
C        new routine PCKEUL.  Otherwise the code calculates
C        the Euler angles from the P_constants kernel file.
C
C-     SPICELIB Version 2.0.0, 04-SEP-1991 (NJB)
C
C         Updated to handle P_constants referenced to different epochs
C         and inertial reference frames.
C
C         The header was updated to specify that the inertial reference
C         frame used by BODMAT is restricted to be J2000.
C
C         BODMAT now checks the kernel pool for presence of the
C         variables
C
C            BODY#_CONSTANTS_REF_FRAME
C
C         and
C
C            BODY#_CONSTANTS_JED_EPOCH
C
C         where # is the NAIF integer code of the barycenter of a
C         planetary system or of a body other than a planet or
C         satellite.  If either or both of these variables are
C         present, the P_constants for BODY are presumed to be
C         referenced to the specified inertial frame or epoch.
C         If the epoch of the constants is not J2000, the input
C         time ET is converted to seconds past the reference epoch.
C         If the frame of the constants is not J2000, the rotation from
C         the P_constants' frame to body-fixed coordinates is
C         transformed to the rotation from J2000 coordinates to
C         body-fixed coordinates.
C
C         For efficiency reasons, this routine now duplicates much
C         of the code of BODEUL so that it doesn't have to call BODEUL.
C         In some cases, BODEUL must covert Euler angles to a matrix,
C         rotate the matrix, and convert the result back to Euler
C         angles.  If this routine called BODEUL, then in such cases
C         this routine would convert the transformed angles back to
C         a matrix.  That would be a bit much....
C
C
C-    Beta Version 1.1.0, 16-FEB-1989 (IMU) (NJB)
C
C        Examples section completed.  Declaration of unused variable
C        FOUND removed.
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
      DOUBLE PRECISION      TSIPM   ( 6, 6 )

      INTEGER               I
      INTEGER               J
 
C
C     Standard SPICE Error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'BODMAT' )
      END IF

C
C     Get 6x6 state transformation from TISBOD. If succeeded, pull out
C     left-top 3x3 matrix.
C
      CALL TISBOD ( 'J2000', BODY, ET, TSIPM )
 
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'BODMAT' )
         RETURN
      END IF

      DO  I = 1, 3
         DO   J  = 1, 3
            TIPM  ( I, J ) = TSIPM  ( I, J )
         END DO
      END DO

      CALL CHKOUT ( 'BODMAT' )
      RETURN
      END
