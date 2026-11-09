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
  $Id: tpih1.c,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
*=============================================================================
*
* tpih1 tests wcspih(), the WCS FITS parser for image headers, and wcsfix(),
* which translates non-standard constructs.  It reads a test header and uses
* wcsprt() to print the resulting wcsprm structs.
*
* Input comes from file "pih.fits" using either fits_hdr2str() from CFITSIO
* if the DO_CFITSIO preprocessor is defined, or read directly using fgets()
* otherwise.
*
*---------------------------------------------------------------------------*/

// Must be included before anything else.
#include <wcsconfig_tests.h>

#include <stdio.h>
#include <string.h>

#if defined HAVE_CFITSIO && defined DO_CFITSIO
#include <fitsio.h>
#endif

#include <wcs.h>
#include <wcshdr.h>
#include <wcsfix.h>
#include <wcsprintf.h>

int main()

{
  // Set line buffering in case stdout is redirected to a file, otherwise
  // stdout and stderr messages will be jumbled (stderr is unbuffered).
  setvbuf(stdout, NULL, _IOLBF, 0);

  printf("Testing WCSLIB parser for FITS image headers (tpih1.c)\n"
         "------------------------------------------------------\n\n");

  // Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
  const char infile[] = "pih.fits";
  int nkeyrec = 0, status = 0;

#if defined HAVE_CFITSIO && defined DO_CFITSIO
  fitsfile *fptr;
  if (fits_open_file(&fptr, infile, READONLY, &status)) {
    fits_report_error(stderr, status);
    return 1;
  }

  char *header;
  if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
    fits_report_error(stderr, status);
    return 1;
  }

  fits_close_file(fptr, &status);
#else
  FILE *fptr;
  if ((fptr = fopen(infile, "r")) == 0) {
    fprintf(stderr, "ERROR opening %s\n", infile);
    return 1;
  }

  char header[288001];
  int hp = 0;
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

      memcpy(header+hp, keyrec, 80);
      hp += 80;
      nkeyrec++;

      if (strncmp(keyrec, "END     ", 8) == 0) {
        // An END keyrecord was read, but read the rest of the block.
        gotend = 1;
      }
    }

    if (gotend) break;
  }
  fclose(fptr);
#endif

  fprintf(stderr, "Found %d non-comment header keyrecords.\n\n", nkeyrec);


  // Parse the header, allowing all recognized non-standard WCS keywords and
  // usage.  All WCS keyrecords are culled from the header, illegal ones are
  // reported.
  fprintf(stderr, "\nIllegal or extraneous WCS header keyrecords rejected "
                  "by wcspih():\n");

  int relax = WCSHDR_all;
  int ctrl  = -2;
  int nreject, nwcs;
  struct wcsprm *wcs;
  if ((status = wcspih(header, nkeyrec, relax, ctrl, &nreject, &nwcs,
                       &wcs))) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcs_errmsg[status]);
  }
  if (!nreject) fprintf(stderr, "(none)\n");


  // List the remaining keyrecords.
  printf("\n\nNon-WCS header keyrecords ignored by wcspih():\n");
  char *hptr = header;
  while (*hptr) {
    printf("%.80s\n", hptr);
    hptr += 80;
  }
#if defined HAVE_CFITSIO && defined DO_CFITSIO
  fits_free_memory(header, &status);
#endif


  // Summarize what was found.
  int alts[27];
  status = wcsidx(nwcs, &wcs, alts);
  printf("\n\nFound %d alternate coordinate descriptions with indices:\n  ",
         nwcs);
  for (int a = 'A'; a <= 'Z'; a++) {
    printf("%2c", a);
  }
  printf("\n");
  for (int ialt = 0; ialt < 27; ialt++) {
    if (alts[ialt] < 0) {
      printf(" -");
    } else {
      printf("%2d", alts[ialt]);
    }
  }
  printf("\n");


  // Fix non-standard usage and print each of the wcsprm structs.  The output
  // from wcsprt() will be written to an internal buffer and then printed just
  // to show that it can be done.
  wcsprintf_set(0x0);
  for (int iwcs = 0; iwcs < nwcs; iwcs++) {
    wcsprintf("\n------------------------------------"
              "------------------------------------\n");

    // Fix non-standard WCS keyvalues.
    int stat[NWCSFIX];
    if ((status = wcsfix(7, 0, wcs+iwcs, stat))) {
      printf("wcsfix ERROR, status returns: (");
      for (int ifix = 0; ifix < NWCSFIX; ifix++) {
        printf(ifix ? ", %d" : "%d", stat[ifix]);
      }
      printf(")\n\n");
    }

    if ((status = wcsset(wcs+iwcs))) {
      fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
      continue;
    }

    if ((status = wcsprt(wcs+iwcs))) {
      fprintf(stderr, "wcsprt ERROR %d: %s.\n", status, wcs_errmsg[status]);
    }
  }
  printf("%s", wcsprintf_buf());

  // Defeat spurious reporting of memory leaks.
  wcsprintf_set(stdout);
  wcsvfree(&nwcs, &wcs);

  return 0;
}
