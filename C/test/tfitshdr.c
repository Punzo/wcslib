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
  $Id: tfitshdr.c,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
*=============================================================================
*
* tfitshdr tests fitshdr(), the FITS parser for image headers, by reading a
* test header and printing the resulting fitskey structs.
*
* Input comes from file 'fitshdr.fits' using either fits_hdr2str() from
* CFITSIO if the DO_CFITSIO preprocessor is defined, or read directly using
* fgets() otherwise.
*
* wcshdr() is called first to extract all WCS-related keyrecords from the
* input header before passing it on to fitshdr().
*
*---------------------------------------------------------------------------*/

// Must be included before anything else.
#include <wcsconfig_tests.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined HAVE_CFITSIO && defined DO_CFITSIO
#include <fitsio.h>
#endif

#include <fitshdr.h>
#include <wcshdr.h>
#include <wcsutil.h>

int main()

{
  // Set line buffering in case stdout is redirected to a file, otherwise
  // stdout and stderr messages will be jumbled (stderr is unbuffered).
  setvbuf(stdout, NULL, _IOLBF, 0);

  printf("Testing FITS image header parser (tfitshdr.c)\n"
         "---------------------------------------------\n\n");

  // Read in the FITS header.
  char infile[] = "fitshdr.fits";
  int  nkeyrec, status = 0;

#if defined HAVE_CFITSIO && defined DO_CFITSIO
  fitsfile *fptr;
  if (fits_open_file(&fptr, infile, READONLY, &status)) {
    fits_report_error(stderr, status);
    return 1;
  }

  char *header;
  if (fits_hdr2str(fptr, 0, NULL, 0, &header, &nkeyrec, &status)) {
    fits_report_error(stderr, status);
    return 1;
  }

  fits_close_file(fptr, &status);
#else
  FILE *fptr;
  if ((fptr = fopen(infile, "r")) == 0) {
    printf("ERROR opening %s\n", infile);
    return 1;
  }

  char keyrec[81], header[288001];
  int k = 0, gotend = 0;
  nkeyrec = 0;
  for (int j = 0; j < 100; j++) {
    for (int i = 0; i < 36; i++) {
      if (fgets(keyrec, 81, fptr) == 0) {
        break;
      }

      memcpy(header+k, keyrec, 80);
      k += 80;
      nkeyrec++;

      if (strncmp(keyrec, "END       ", 10) == 0) {
        // An END keyrecord was read, but read the rest of the block.
        gotend = 1;
      }
    }

    if (gotend) break;
  }
  fclose(fptr);
#endif

  printf("Found %d header keyrecords.\n", nkeyrec);


  // Cull recognized, syntactically valid WCS keyrecords from the header.
  int relax = WCSHDR_all;
  int ctrl = -1;
  int nreject, nwcs;
  struct wcsprm *wcs;
  if ((status = wcspih(header, nkeyrec, relax, ctrl, &nreject, &nwcs,
                       &wcs))) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcs_errmsg[status]);
    return 1;
  }
  wcsvfree(&nwcs, &wcs);

  // Number remaining.
  nkeyrec = strlen(header) / 80;


  // Specific keywords to be located or culled.
  int nkeyids = 7;
  struct fitskeyid keyids[8];
  strncpy(keyids[0].name, "SIMPLE  ", 12);
  strncpy(keyids[1].name, "BITPIX  ", 12);
  strncpy(keyids[2].name, "NAXIS   ", 12);
  strncpy(keyids[3].name, "COMMENT ", 12);
  strncpy(keyids[4].name, "HISTORY ", 12);
  strncpy(keyids[5].name, "        ", 12);
  strncpy(keyids[6].name, "END     ", 12);

  if (nkeyids) {
    printf("\nThe following keyrecords will not be listed:\n");
    for (int i = 0; i < nkeyids; i++) {
      printf("  \"%8s\"\n", keyids[i].name);
    }
  }


  // Parse the header.
  struct fitskey *keys;
  if ((status = fitshdr(header, nkeyrec, nkeyids, keyids, &nreject, &keys))) {
    printf("fitskey ERROR %d: %s.\n", status, fitshdr_errmsg[status]);
  }
#if defined HAVE_CFITSIO && defined DO_CFITSIO
  fits_free_memory(header, &status);
#endif

  // Report the results.
  printf("\n%d header keyrecords parsed by fitshdr(), %d rejected:\n\n",
    nkeyrec, nreject);
  struct fitskey *kptr = keys;
  for (int i = 0; i < nkeyrec; i++, kptr++) {
    // Skip syntactically valid keyrecords that were indexed.
    if (kptr->keyno < 0 && !kptr->status) continue;

    // Basic keyrecord info.
    printf("%4d%5d  %-8s%3d", kptr->keyno, kptr->status, kptr->keyword,
                              kptr->type);

    // Format the keyvalue for output.
    char text[80];
    switch (abs(kptr->type)%10) {
    case 1:
      // Logical.
      snprintf(text, 80, "%c", kptr->keyvalue.i?'T':'F');
      break;
    case 2:
      // 32-bit signed integer.
      snprintf(text, 80, "%d", kptr->keyvalue.i);
      break;
    case 3:
      // 64-bit signed integer.
#ifdef WCSLIB_INT64
         snprintf(text, 80, "%+lld", kptr->keyvalue.k);
#else
         if (kptr->keyvalue.k[2]) {
           snprintf(text, 80, "%+d%09d%09d", kptr->keyvalue.k[2],
                                         abs(kptr->keyvalue.k[1]),
                                         abs(kptr->keyvalue.k[0]));
         } else {
           snprintf(text, 80, "%+d%09d",     kptr->keyvalue.k[1],
                                         abs(kptr->keyvalue.k[0]));
         }
#endif
       break;
    case 4:;
      // Very long integer.
      int k = 0;
      for (int j = 7; j > 0; j--) {
        if (kptr->keyvalue.l[j]) {
          k = j;
          break;
        }
      }

      snprintf(text, 80, "%+d", kptr->keyvalue.l[k]);
      char *tp = text + strlen(text);
      for (int j = k-1; j >= 0; j--) {
        size_t tsize = (text + 80) - tp;
        snprintf(tp, tsize, "%09d", abs(kptr->keyvalue.l[j]));
        tp += strlen(tp);
      }

      break;
    case 5:
      // Float.
      snprintf(text, 80, "%+13.6e", kptr->keyvalue.f);
      break;
    case 6:
      // Int complex.
      snprintf(text, 80, "%.0f  %.0f", kptr->keyvalue.c[0],
                                       kptr->keyvalue.c[1]);
      break;
    case 7:
      // Float complex.
      snprintf(text, 80, "%+13.6e  %+13.6e", kptr->keyvalue.c[0],
                                             kptr->keyvalue.c[1]);
      break;
    case 8:
      // String.
      snprintf(text, 80, "\"%s\"", kptr->keyvalue.s);
      break;
    default:
      // No value.
      *text = '\0';
      break;
    }

    if (kptr->type > 0) {
      // Keyvalue successfully extracted.
      printf("  %s", text);
    } else if (kptr->type < 0) {
      // Syntax error of some type while extracting the keyvalue.
      printf("  (%s)", text);
    }

    // Units?
    if (kptr->ulen) {
      printf(" %.*s", kptr->ulen-2, kptr->comment+1);
    }

    // Comment text or reject keyrecord.
    printf("\n%s\n", kptr->comment);
  }


  // Print indexes.
  printf("\n\nIndexes of selected keywords:\n");
  for (int i = 0; i < nkeyids; i++) {
    printf("%-8s%5d%5d%5d", keyids[i].name, keyids[i].count, keyids[i].idx[0],
      keyids[i].idx[1]);

    // Print logical (SIMPLE) and integer (BITPIX, NAXIS) values.
    if (keyids[i].count) {
      kptr = keys + keyids[i].idx[0];
      printf("%4d", kptr->type);

      if (kptr->type == 1) {
        printf("%5c", kptr->keyvalue.i?'T':'F');
      } else if (kptr->type == 2) {
        printf("%5d", kptr->keyvalue.i);
      }
    }
    printf("\n");
  }

  wcsdealloc(keys);

  return 0;
}
