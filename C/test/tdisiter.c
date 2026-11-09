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
  $Id: tdisiter.c,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
*=============================================================================
*
* tdisiter tests the external setting of ITERMAX in disx2p(), particularly
* for distortions such as SIP that may optionally provide an approximation
* of their inverses.
*
*---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <wcserr.h>
#include <wcshdr.h>
#include <wcsprintf.h>
#include <dis.h>

int disitermax(int itermax);

#define HEADER_SIZE 36000


int main(int argc, char *argv[])

{
  // Set line buffering in case stdout is redirected to a file, otherwise
  // stdout and stderr messages will be jumbled (stderr is unbuffered).
  setvbuf(stdout, NULL, _IOLBF, 0);

  wcserr_enable(1);
  wcsprintf_set(stdout);

  // Optional file name specified?
  const char *infile = "SIP.fits";
  if (1 < argc) {
    infile = argv[1];
  }


  // Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
  FILE *fptr;
  if ((fptr = fopen(infile, "r")) == 0) {
    wcsprintf("ERROR opening %s\n", infile);
    return 1;
  }

  int naxis[4];
  memset(naxis, 0, 2*sizeof(int));

  char header[288001];
  int hp = 0;
  int nkeyrec = 0;
  int gotend = 0;
  for (int iblock = 0; iblock < 100; iblock++) {
    for (int ikeyrec = 0; ikeyrec < 36; ikeyrec++) {
      char keyrec[81];
      if (fgets(keyrec, 81, fptr) == 0) {
        break;
      }

      if (strncmp(keyrec, "        ", 8) == 0) continue;
      if (strncmp(keyrec, "COMMENT ", 8) == 0) continue;
      if (strncmp(keyrec, "HISTORY ", 8) == 0) continue;

      if (strncmp(keyrec, "NAXIS", 5) == 0) {
        if (keyrec[5] == ' ') {
          int n;
          sscanf(keyrec+10, "%d", &n);
          if (4 < n) {
            wcsprintf("ERROR, can't handle more than 4 axes.\n");
            return 1;
          }
          continue;
        }

        int j, n;
        sscanf(keyrec+5, "%d = %d", &j, &n);
        naxis[j-1] = n;
        continue;
      }

      strncpy(header+hp, keyrec, 80);
      hp += 80;
      nkeyrec++;

      if (strncmp(keyrec, "END       ", 10) == 0) {
        // An END keyrecord was read, but read the rest of the block.
        gotend = 1;
      }
    }

    if (gotend) break;
  }
  fclose(fptr);


  // Parse the header.
  int nreject, nwcs;
  struct wcsprm *wcs;
  if ((wcspih(header, nkeyrec, WCSHDR_none, 2, &nreject, &nwcs, &wcs))) {
    wcsperr(wcs, 0x0);
    return 1;
  }

  if (wcsset(wcs)) {
    wcsperr(wcs, 0x0);
    return 1;
  }

  struct disprm *dis = wcs->lin.dispre;

  disprt(dis);
  wcsprintf("\n");

  double rawcrd[] = {1.0, 1.0, 1.0, 1.0};
  double discrd[4];
  if (disp2x(dis, rawcrd, discrd)) {
    disperr(dis, 0x0);
    return 1;
  }

  wcsprintf("    pix: %18.12f %18.12f\n", rawcrd[0], rawcrd[1]);
  wcsprintf(" -> img: %18.12f %18.12f\n", discrd[0], discrd[1]);

  disitermax(32);
  if (disx2p(dis, discrd, rawcrd)) {
    disperr(dis, 0x0);
    return 1;
  }

  wcsprintf(" -> pix: %18.12f %18.12f  (ITERMAX == %d)\n",
    rawcrd[0], rawcrd[1], disitermax(-1));

  disitermax(0);
  if (disx2p(dis, discrd, rawcrd)) {
    disperr(dis, 0x0);
    return 1;
  }

  wcsprintf(" -> pix: %18.12f %18.12f  (ITERMAX == %d)\n",
    rawcrd[0], rawcrd[1], disitermax(-1));


  wcsvfree(&nwcs, &wcs);

  return 0;
}
