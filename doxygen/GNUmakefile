#-----------------------------------------------------------------------------
# GNU makefile for building the documentation for WCSLIB 8.10
#
# Summary of the main targets
# ---------------------------
#   build:     Invoke 'doxygen' etc. to compile the documentation.
#   clean:     Delete rubbish files.
#   cleaner (distclean or realclean): clean, and also delete intermediate
#              files.
#   cleanest:  cleaner, and also delete the generated documentation.
#
# Notes:
#   1) If you need to make changes then preferably modify ../makedefs.in
#      instead and re-run configure.
#
# Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
# http://www.atnf.csiro.au/computing/software/wcs
# $Id: GNUmakefile,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
#-----------------------------------------------------------------------------
# Get configure settings.
SUBDIR := doxygen
include ../makedefs

DOXHDRS  := $(notdir $(wildcard ../C/*.h))
DOXMARX  := $(DOXHDRS:.h=.dox)
EXTRADOX := mainpage.dox $(wildcard *_extras.dox)

# Pattern rules
#--------------

%.dox : ../C/%.h doxextr doxextr.sed %.sed Doxyfile
	-@ echo ''
	 @ if [ -f $*.sed ] ; then \
	     echo "doxextr $< | sed -f $*.sed | sed -f doxextr.sed > $@" ; \
	     ./doxextr $< | sed -f $*.sed | sed -f doxextr.sed > $@ ; \
	   else \
	     echo "doxextr $< | sed -f doxextr.sed > $@" ; \
	     ./doxextr $< | sed -f doxextr.sed > $@ ; \
	   fi

%.c : %.l
	-@ echo ''
	-@ $(RM) $@
	   $(FLEX) $(FLFLAGS) -t $< | sed -e 's/^[	 ]*#/#/' > $@

% : %.c
	-@ echo ''
	   $(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $<

# Don't fret if there's no specific sed script.
%.sed : ;

# Static rules
#-------------

.PHONY : build clean cleaner cleanest distclean html install pdf realclean \
         uninstall

build : doxextr html pdf

html : ../html/index.html

pdf : ../wcslib.pdf

../wcslib.pdf : latex/refman.tex
	-  $(MAKE) -C latex refman.pdf
	-@ $(RM) $@
	   mv latex/refman.pdf $@

../html/index.html latex/refman.tex : doxextr wcslib.h
	-@ echo ''
	   doxygen
	   sed -e 's/let~}}/let~}\\relax}/' latex/doxygen.sty > doxygen.sty
	   mv latex/doxygen.sty latex/doxygen.sty.bak
	   mv doxygen.sty latex/

wcslib.h : $(DOXMARX) $(EXTRADOX)
	 @ cat ../C/wcslib.h *.dox > $@

install :
	-  if [ ! -d "$(PDFDIR)" ] ; then \
	     $(INSTALL) -d -m 775 $(PDFDIR) ; \
	   fi
	   $(INSTALL) -m 444 ../wcslib.pdf $(PDFDIR)
	-  if [ ! -d "$(HTMLDIR)/html/search" ] ; then \
	     $(INSTALL) -d -m 775 $(HTMLDIR)/html/search ; \
	   fi
	 @ mv ../html/search ../html/.search
	   $(INSTALL) -m 444 ../html/* $(HTMLDIR)/html
	 @ mv ../html/.search ../html/search
	   $(INSTALL) -m 444 ../html/search/* $(HTMLDIR)/html/search

uninstall :
	-  $(RM) $(PDFDIR)
	-  $(RM) $(HTMLDIR)

clean :
	- $(RM) core

cleaner distclean realclean : clean
	- $(RM) *.log
	- $(RM) doxextr doxextr.c wcslib.h
	- $(RM) $(DOXMARX)
	- $(RM) wcslib.tag
	- $(RM) -r ./latex

cleanest : cleaner
	- $(RM) -r ../html ../wcslib.pdf

GNUmakefile : ../makedefs ;

../makedefs : ../makedefs.in
	   cd .. && ./config.status
