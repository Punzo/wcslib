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
  $Id: WCSlogo.c,v 8.10 2026/09/27 05:44:35 mcalabre Exp $
*=============================================================================
*
* Draws the WCSLIB logo.  Derived from tcel1.
*
*---------------------------------------------------------------------------*/

#include <cpgplot.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cel.h>

int main()

{
  // Initialize.
  struct celprm native;
  celini(&native);

  // Reference angles for the native graticule (in fact, the defaults).
  native.ref[0] = 0.0;
  native.ref[1] = 0.0;

  // Set up Bonne's projection with conformal latitude at +35.
  strncpy(native.prj.code, "BON", 4);
  native.prj.pv[1] = 35.0;

  // Other parameters for the logo.
  const int crval1  =  0;
  const int crval2  = 45;
  const int lonpole = 45;
  const int latpole =  0;

  // Logo size in mm (approximate).
  const float size = 100.0f;


  // Celestial graticule.
  struct celprm celestial;
  celini(&celestial);
  celestial.prj = native.prj;

  celestial.ref[0] = (double)crval1;
  celestial.ref[1] = (double)crval2;
  celestial.ref[2] = (double)lonpole;
  celestial.ref[3] = (double)latpole;


  // PGPLOT initialization.
  char devtyp[80];
  strncpy(devtyp, "WCSlogo.eps/vcps", 80);
  cpgbeg(0, devtyp, 1, 1);
  cpgpap(size/25.4f, 1.0f);

  // Define PGPLOT viewport.
  cpgpage();
  cpgsvp(0.0f, 1.0f, 0.0f, 1.0f);
  cpgswin(-135.0f, 135.0f, -125.0f, 145.0f);

  // Define colours.
  cpgscr(0, 0.2f, 0.2f, 0.4f);  // Background colour.
  cpgscr(1, 1.0f, 1.0f, 0.0f);
  cpgscr(2, 1.0f, 1.0f, 1.0f);
  cpgscr(3, 0.5f, 0.5f, 0.8f);
  cpgscr(4, 0.8f, 0.5f, 0.5f);
  cpgscr(5, 0.8f, 0.8f, 0.8f);
  cpgscr(6, 0.5f, 0.5f, 0.8f);
  cpgscr(7, 0.8f, 0.5f, 0.5f);
  cpgscr(8, 0.3f, 0.5f, 0.3f);  // Native grid colour.


  int    stat[361];
  double lat[181], lng[361], phi[361], theta[361], x[361], y[361];

  // Do black background.
  lng[0] = -179.99;
  for (int j = 0, ilat = -90; ilat <= 90; ilat++, j++) {
    lat[j] = (double)ilat;
  }

  cels2x(&native, 1, 181, 1, 1, lng, lat, phi, theta, x, y, stat);

  lng[0] = 179.99;
  for (int j = 0, ilat = 90; ilat >= -90; ilat--, j++) {
    lat[j] = (double)ilat;
  }

  cels2x(&native, 1, 181, 1, 1, lng, lat, phi, theta, x+180, y+180, stat);

  float xr[512], yr[512];
  for (int j = 0; j <= 360; j++) {
    xr[j] = -1.03f * (float)x[j];
    yr[j] =  1.03f * (float)y[j];
  }

  cpgsci(0);
  cpgpoly(361, xr, yr);


  // Draw the native graticule faintly in the background.
  cpgsci(8);

  // Draw native meridians of longitude.
  for (int j = 0, ilat = -90; ilat <= 90; ilat++, j++) {
    lat[j] = (double)ilat;
  }

  for (int ilng = -180; ilng <= 180; ilng += 15) {
    lng[0] = (double)ilng;
    if (ilng == -180) lng[0] = -179.99;
    if (ilng ==  180) lng[0] =  179.99;

    // Dash the longitude of the celestial pole.
    if ((ilng-lonpole)%360 == 0) {
      cpgsls(2);
      cpgslw(5);
    }

    cels2x(&native, 1, 181, 1, 1, lng, lat, phi, theta, x, y, stat);

    int k = 0;
    for (int j = 0; j < 181; j++) {
      if (stat[j]) {
        if (k > 1) cpgline(k, xr, yr);
        k = 0;
        continue;
      }

      xr[k] = -(float)x[j];
      yr[k] =  (float)y[j];
      k++;
    }

    cpgline(k, xr, yr);
    cpgsls(1);
    cpgslw(1);
  }

  // Draw native parallels of latitude.
  lng[0]   = -179.99;
  lng[360] =  179.99;
  for (int j = 1, ilng = -179; ilng < 180; ilng++, j++) {
    lng[j] = (double)ilng;
  }

  for (int ilat = -90; ilat <= 90; ilat += 15) {
    lat[0] = (double)ilat;

    cels2x(&native, 361, 1, 1, 1, lng, lat, phi, theta, x, y, stat);

    int k = 0;
    for (int j = 0; j < 361; j++) {
      if (stat[j]) {
        if (k > 1) cpgline(k, xr, yr);
        k = 0;
        continue;
      }

      xr[k] = -(float)x[j];
      yr[k] =  (float)y[j];
      k++;
    }

    cpgline(k, xr, yr);
  }


  // Draw a colour-coded celestial coordinate graticule.
  int ci = 1;

  // Draw celestial meridians of longitude.
  for (int j = 0, ilat = -90; ilat <= 90; ilat++, j++) {
    lat[j] = (double)ilat;
  }

  for (int ilng = -180; ilng <= 180; ilng += 15) {
    lng[0] = (double)ilng;

    if (++ci > 7) ci = 2;
    cpgsci(ilng?ci:1);

    // Dash the reference longitude.
    if ((ilng-crval1)%360 == 0) {
      cpgsls(2);
      cpgslw(5);
    }

    cels2x(&celestial, 1, 181, 1, 1, lng, lat, phi, theta, x, y, stat);

    int k = 0;
    for (int j = 0; j < 181; j++) {
      if (stat[j]) {
        if (k > 1) cpgline(k, xr, yr);
        k = 0;
        continue;
      }

      // Test for discontinuities.
      if (j > 0) {
        if (fabs(x[j]-x[j-1]) > 4.0 || fabs(y[j]-y[j-1]) > 4.0) {
          if (k > 1) cpgline(k, xr, yr);
          k = 0;
        }
      }

      xr[k] = -(float)x[j];
      yr[k] =  (float)y[j];
      k++;
    }

    cpgline(k, xr, yr);
    cpgsls(1);
    cpgslw(1);
  }

  // Draw celestial parallels of latitude.
  for (int j = 0, ilng = -180; ilng <= 180; ilng++, j++) {
    lng[j] = (double)ilng;
  }

  ci = 1;
  for (int ilat = -90; ilat <= 90; ilat += 15) {
    lat[0] = (double)ilat;

    if (++ci > 7) ci = 2;
    cpgsci(ilat?ci:1);

    // Dash the reference latitude.
    if (ilat == crval2) {
      cpgsls(2);
      cpgslw(5);
    }

    cels2x(&celestial, 361, 1, 1, 1, lng, lat, phi, theta, x, y, stat);

    int k = 0;
    for (int j = 0; j < 361; j++) {
      if (stat[j]) {
        if (k > 1) cpgline(k, xr, yr);
        k = 0;
        continue;
      }

      // Test for discontinuities.
      if (j > 0) {
        if (fabs(x[j]-x[j-1]) > 4.0 || fabs(y[j]-y[j-1]) > 4.0) {
          if (k > 1) cpgline(k, xr, yr);
          k = 0;
        }
      }

      xr[k] = -(float)x[j];
      yr[k] =  (float)y[j];
      k++;
    }

    cpgline(k, xr, yr);
    cpgsls(1);
    cpgslw(1);
  }

  cpgend();

  printf("WCSLIB logo written to WCSlogo.eps.\n");

  return 0;
}
