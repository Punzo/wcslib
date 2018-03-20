/* wcsconfig_f77.h.  Generated from wcsconfig_f77.h.in by configure.  */
/*============================================================================
*
* wcsconfig_f77.h is generated from wcsconfig_f77.h.in by 'configure'.  It
* contains C preprocessor definitions for building the WCSLIB 5.18 Fortran
* wrappers.
*
* Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
* http://www.atnf.csiro.au/people/Mark.Calabretta
* $Id: wcsconfig_f77.h.in,v 5.18 2018/01/10 08:32:14 mcalabre Exp $
*===========================================================================*/

/* Integer array type large enough to hold an address.  Set here to int[2] for
 * 64-bit addresses, but could be defined as int* on 32-bit machines. */
typedef int iptr[2];

/* Macro for mangling Fortran subroutine names that do not contain
 * underscores.  Typically a name like "WCSINI" (case-insensitive) will become
 * something like "wcsini_" (case-sensitive).  The Fortran wrappers, which are
 * written in C, are preprocessed into names that match the latter.  The macro
 * takes two arguments which specify the name in lower and upper case. */
#define F77_FUNC(name,NAME) name ## _
