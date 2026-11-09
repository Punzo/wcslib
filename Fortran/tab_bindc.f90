!=============================================================================
! WCSLIB 8.10 - an implementation of the FITS WCS standard.
! Copyright (C) 1995-2026, Mark Calabretta
!
! This file is part of WCSLIB.
!
! WCSLIB is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with WCSLIB.  If not, see http://www.gnu.org/licenses.
!
! Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
! http://www.atnf.csiro.au/computing/software/wcs
! $Id: tab_bindc.f90,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
!=============================================================================

INTEGER FUNCTION TABPERR (TAB, PREFIX)
  INTEGER :: TAB(*)
  CHARACTER :: PREFIX(*)

  INTERFACE
    INTEGER(C_INT) FUNCTION TABPERR_C (TAB, PREFIX) BIND (C)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(C_INT) :: TAB(*)
      CHARACTER(KIND=C_CHAR, LEN=1) :: PREFIX(72)
    END FUNCTION TABPERR_C
  END INTERFACE

  TABPERR = TABPERR_C (TAB, PREFIX)
END FUNCTION TABPERR
