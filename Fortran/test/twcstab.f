*=======================================================================
*
* WCSLIB 5.15 - an implementation of the FITS WCS standard.
* Copyright (C) 1995-2016, Mark Calabretta
*
* This file is part of WCSLIB.
*
* WCSLIB is free software: you can redistribute it and/or modify it
* under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* WCSLIB is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
* License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with WCSLIB.  If not, see http://www.gnu.org/licenses.
*
* Direct correspondence concerning WCSLIB to mark@calabretta.id.au
*
* Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
* http://www.atnf.csiro.au/people/Mark.Calabretta
* $Id: twcstab.f,v 5.15 2016/04/05 12:55:12 mcalabre Exp $
*=======================================================================

      PROGRAM TWCSTAB
*-----------------------------------------------------------------------
*
* TWCSTAB tests WCSTAB and also provides sample code for using it in
* conjunction with WCSPIH and FTWCST.  Although this example and FTWCST
* are based on the CFITSIO library, WCSTAB itself is completely
* independent of it.
*
* We assume that the input file, ../C/wcstab.fits, has already been
* generated by running the C version of twcstab.
*
* WCSP and WTB, which are meant to hold addresses, are declared as
* INTEGER arrays of length 2 to accomodate 64-bit machines for which
* sizeof(void *) = 2*sizeof(int).
*=======================================================================

      LOGICAL   GOTEND
      INTEGER   BLOKSZ, I, J, K, IERR, IUNIT, NKEYRC, NREJECT, NWCS,
     :          NWTB, STATUS, WCSP(2), WTB(2)
      CHARACTER KEYREC*80, HEADER*28801, INFILE*16

*     On some systems, such as Sun Sparc, the struct MUST be aligned
*     on a double precision boundary, done here using an equivalence.
*     Failure to do this may result in mysterious "bus errors".
      INCLUDE 'wcs.inc'
      INCLUDE 'wcshdr.inc'
      INCLUDE 'wcsfix.inc'
      INCLUDE 'getwcstab.inc'
      INTEGER STAT(WCSFIX_NWCS)
      INTEGER WCS(WCSLEN)
      DOUBLE PRECISION DUMMY
      EQUIVALENCE (WCS,DUMMY)

      DATA INFILE /'../C/wcstab.fits'/
*-----------------------------------------------------------------------
      WRITE (*, 10)
 10   FORMAT ('Testing WCSTAB and associated routines (twcstab.f)',/,
     :        '--------------------------------------------------',/)

*     Open the FITS test file.
      IUNIT = 1
      STATUS = 0
      CALL FTOPEN (IUNIT, INFILE, 0, BLOKSZ, STATUS)
      IF (STATUS.NE.0) THEN
        CALL FLUSH(6)
        CALL FTRPRT ('STDERR', STATUS)
        GO TO 999
      END IF

*     Read the primary header; unfortunately there is no FITSIO
*     equivalent of CFITSIO's fits_hdr2str() so do it the long way.
      OPEN (UNIT=1, FILE=INFILE, FORM='FORMATTED', ACCESS='DIRECT',
     :      RECL=80, IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE (*, 20) IERR, INFILE
 20     FORMAT ('ERROR',I3,' opening ',A)
        GO TO 999
      END IF

*     Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
      K = 1
      NKEYRC = 0
      GOTEND = .FALSE.
      DO 50 J = 0, 100
        DO 40 I = 1, 36
          READ (1, '(A80)', REC=36*J+I, IOSTAT=IERR) KEYREC
          IF (IERR.NE.0) THEN
            WRITE (*, 30) IERR
 30         FORMAT ('ERROR',I3,' reading header.')
            GO TO 999
          END IF

          IF (KEYREC(:8).EQ.'        ') GO TO 40
          IF (KEYREC(:8).EQ.'COMMENT ') GO TO 40
          IF (KEYREC(:8).EQ.'HISTORY ') GO TO 40

          HEADER(K:) = KEYREC
          K = K + 80
          NKEYRC = NKEYRC + 1

          IF (KEYREC(:8).EQ.'END     ') THEN
*           An END keyrecord was read, read the rest of the block.
            GOTEND = .TRUE.
          END IF
 40     CONTINUE

        IF (GOTEND) GO TO 60
 50   CONTINUE

 60   CLOSE (UNIT=1)

*-----------------------------------------------------------------------
* Basic steps required to interpret a FITS WCS header, including -TAB.
*-----------------------------------------------------------------------

*     Parse the primary header of the FITS file.
      STATUS = WCSPIH (HEADER, NKEYRC, WCSHDR_all, 2, NREJECT, NWCS,
     :                 WCSP)
      IF (STATUS.NE.0) THEN
        WRITE (*, 70) STATUS, WCSHDR_ERRMSG(STATUS)
 70     FORMAT ('WCSPIH ERROR',I2,A)
        GO TO 999
      END IF

*     Copy into our WCSPRM struct.
      IERR = WCSVCOPY (WCSP, 0, WCS)

*     Read coordinate arrays from the binary table extension.
      STATUS = WCSGET (WCS, WCS_NWTB, NWTB)
      STATUS = WCSGET (WCS, WCS_WTB,  WTB)
      STATUS = FTWCST (IUNIT, NWTB, WTB, STATUS)
      IF (STATUS.NE.0) THEN
        CALL FLUSH(6)
        CALL FTRPRT ('STDERR', STATUS)
        GO TO 999
      END IF

*     Fix non-standard WCS keyvalues.
      STATUS = WCSFIX (7, 0, WCS, STAT)
      IF (STATUS.NE.0) THEN
        WRITE (*, 80) (STAT(I), I=1,WCSFIX_NWCS)
 80     FORMAT ('WCSFIX ERROR, status returns: ',10(I2,:,','))
        GO TO 999
      END IF

*-----------------------------------------------------------------------
* The wcsprm struct is now ready for use.
*-----------------------------------------------------------------------

*     Finished with the FITS file.
      CALL FTCLOS (IUNIT, STATUS)

*     Initialize the wcsprm struct, taking memory allocated by FTWCST.
      IF (STATUS.EQ.0) STATUS = WCSSET (WCS)
      IF (STATUS.NE.0) THEN
        WRITE (*, 90) STATUS, WCS_ERRMSG(STATUS)
 90     FORMAT ('WCSSET ERROR',I2,A)
        GO TO 998
      END IF

*     Do something with it.
      CALL FLUSH(6)
      STATUS = WCSPRT (WCS)

*     Clean up.
 998  STATUS = WCSFREE (WCS)
      STATUS = WCSVFREE (NWCS, WCSP)

 999  CONTINUE
      END