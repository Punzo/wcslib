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
  $Id: tdis2.c,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
*=============================================================================
*
* tdis2 tests wcssub()'s handling of distortion functions.  Input comes from
* the FITS file specified as an argument, or else from SIP.fits.
*
*---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <wcserr.h>
#include <wcshdr.h>
#include <wcsprintf.h>
#include <lin.h>
#include <wcs.h>

#define HEADER_SIZE 36000

int main(int argc, char *argv[])

{
  // Set line buffering in case stdout is redirected to a file, otherwise
  // stdout and stderr messages will be jumbled (stderr is unbuffered).
  setvbuf(stdout, NULL, _IOLBF, 0);

  wcserr_enable(1);
  wcsprintf_set(stdout);

  wcsprintf("Testing wcssub() with distortions (tdis2.c)\n"
            "-------------------------------------------\n");

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

      memcpy(header+hp, keyrec, 80);
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

  // Extract the coordinate description for a transposed subimage and prepend
  // a new axis.  Also tests wcssub() on a struct that hasn't been set up.
  int nsub = 3;
  int axes[4];
  axes[0] = 0;
  axes[1] = WCSSUB_LATITUDE;
  axes[2] = WCSSUB_LONGITUDE;

  int status;
  struct wcsprm wcsext;
  wcsext.flag = -1;
  if ((status = wcssub(1, wcs, &nsub, axes, &wcsext))) {
    wcsperr(&wcsext, "");
    goto cleanup;
  } else if (nsub == 0) {
    printf("None of the requested subimage axes were found.\n");
    goto cleanup;
  }


  // Print the original and extracted structs.
  printf("\nInitial contents of wcsprm struct:\n");
  if ((status = wcsset(wcs))) {
    wcsperr(wcs, "");
    goto cleanup;
  };

  wcsprt(wcs);

  printf("\n\nExtracted contents of wcsprm struct:\n");
  if ((status = wcsset(&wcsext))) {
    wcsperr(&wcsext, "");
    goto cleanup;
  }

  wcsprt(&wcsext);


  // Compute distortion statistics in the initial struct.
  double pixblc[]  = {  1.0,   1.0, 1.0, 1.0};
  double pixtrc[]  = {256.0, 256.0, 1.0, 1.0};
  double pixsamp[] = {  1.0,   1.0, 1.0, 1.0};

  int nsamp;

  double stats[15];
  double *maxdis = stats;
  double *maxtot = maxdis + 4;
  double *avgdis = maxtot + 1;
  double *avgtot = avgdis + 4;
  double *rmsdis = avgtot + 1;
  double *rmstot = rmsdis + 4;

  struct linprm *lin = &(wcs->lin);
  if (linwarp(lin, pixblc, pixtrc, pixsamp, &nsamp,
              maxdis, maxtot, avgdis, avgtot, rmsdis, rmstot)) {
    linperr(lin, 0x0);
    return 1;
  }

  for (int k = 0; k < 12; k++) {
    if (fabs(stats[k]) < 0.0005) stats[k] = 0.0;
  }

  wcsprintf("\n\n"
    "Initial linwarp() statistics computed over %d sample points:\n"
    "  Max distortion, axis 1: %10.5f pixels\n"
    "                  axis 2: %10.5f pixels\n"
    "                  axis 3: %10.5f pixels\n"
    "                  axis 4: %10.5f pixels\n"
    "                   total: %10.5f pixels\n"
    " Mean distortion, axis 1: %10.5f pixels\n"
    "                  axis 2: %10.5f pixels\n"
    "                  axis 3: %10.5f pixels\n"
    "                  axis 4: %10.5f pixels\n"
    "                   total: %10.5f pixels\n"
    "  RMS distortion, axis 1: %10.5f pixels\n"
    "                  axis 2: %10.5f pixels\n"
    "                  axis 3: %10.5f pixels\n"
    "                  axis 4: %10.5f pixels\n"
    "                   total: %10.5f pixels\n",
    nsamp, maxdis[0], maxdis[1], maxdis[2], maxdis[3], *maxtot,
           avgdis[0], avgdis[1], avgdis[2], avgdis[3], *avgtot,
           rmsdis[0], rmsdis[1], rmsdis[2], rmsdis[3], *rmstot);


  // Compute distortion statistics in the extracted struct.
  pixtrc[0]  =   1.0;
  pixtrc[1]  = 256.0;
  pixtrc[2]  = 256.0;

  lin = &(wcsext.lin);
  if (linwarp(lin, pixblc, pixtrc, pixsamp, &nsamp,
              maxdis, maxtot, avgdis, avgtot, rmsdis, rmstot)) {
    linperr(lin, 0x0);
    return 1;
  }

  for (int k = 0; k < 12; k++) {
    if (fabs(stats[k]) < 0.0005) stats[k] = 0.0;
  }

  wcsprintf("\n"
    "linwarp() statistics for extract computed over the same %d points:\n"
    "  Max distortion, axis 1: %10.5f pixels\n"
    "                  axis 2: %10.5f pixels\n"
    "                  axis 3: %10.5f pixels\n"
    "                   total: %10.5f pixels\n"
    " Mean distortion, axis 1: %10.5f pixels\n"
    "                  axis 2: %10.5f pixels\n"
    "                  axis 3: %10.5f pixels\n"
    "                   total: %10.5f pixels\n"
    "  RMS distortion, axis 1: %10.5f pixels\n"
    "                  axis 2: %10.5f pixels\n"
    "                  axis 3: %10.5f pixels\n"
    "                   total: %10.5f pixels\n",
    nsamp, maxdis[0], maxdis[1], maxdis[2], *maxtot,
           avgdis[0], avgdis[1], avgdis[2], *avgtot,
           rmsdis[0], rmsdis[1], rmsdis[2], *rmstot);


cleanup:
  wcsvfree(&nwcs, &wcs);
  wcsfree(&wcsext);

  return status;
}
