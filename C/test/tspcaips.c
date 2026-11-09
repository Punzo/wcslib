/*============================================================================
  WCSLIB 8.10 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2026, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/computing/software/wcs
  $Id: tspcaips.c,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
*=============================================================================
*
* tspcaips does a quick test of spcaips().  Not part of the official test
* suite.
*
*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <spc.h>

int main()

{
  const char *(ctypes[]) = {"FREQ", "VELO", "FELO"};
  const char *(frames[]) = {"-LSR", "-HEL", "-OBS", "    "};

  char ctype[9], ctypeA[9], specsys[9];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      snprintf(ctypeA, 9, "%s%s", ctypes[i], frames[j]);

      for (int v1 = 0; v1 <= 8; v1++) {
        int velref = v1;
        for (int v2 = 0; v2 < 3; v2++) {
          int status = spcaips(ctypeA, velref, ctype, specsys);
          printf("'%s'  %3d  %2d  '%s'  '%s'\n", ctypeA, velref,
            status, ctype, specsys);
          velref += 256;
        }
      }

      printf("\n");
    }
  }

  return 0;
}
