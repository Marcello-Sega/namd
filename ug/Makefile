#####
# define version
#####
NAMD_VERSION = 2.4

#####
# specific names of programs used in make commands
#####
CCPP		= CC
CC		= cc
CD		= cd
COPY		= cp
MOVE		= mv -f
DELETE		= rm -f
ECHO		= echo
TAR		= tar
COMPRESS	= compress
AR		= ar rsc
RANLIB		= /bin/true
COCMD		= co
SHELL 		= /bin/sh


DOCDIR		= .
MAIN_UG_FILE	= $(DOCDIR)/ug

SRCFILES = namd_begindoc.tex namd_copyright.tex namd_macros.tex \
namd_title.tex namd_version.tex ug.tex ug_add.tex ug_avail.tex ug_basic.tex \
ug_intro.tex ug_io.tex ug_macros.tex ug_runit.tex ug_sample.tex ug_start.tex \
ug_xplor.tex ug_fenergy.tex ug_psfgen.tex psfgen_macros.tex ug_alchemy.tex

PDFFIGS = figures/fmaOn.pdf figures/pairlistdist.pdf \
figures/shifting.pdf figures/switching.pdf figures/dual_top.pdf

### latex the documentation
doc:	$(DOCDIR)/ug.ps 

### run latex2html
html:	$(DOCDIR)/ug.ps
	latex2html -local_icons -split 5 -toc_depth 4 -toc_stars \
	-t "NAMD $(NAMD_VERSION) User's Guide" -address "namd@ks.uiuc.edu" \
	-up_url "../../" -up_title "Main NAMD Page" ug

### run pdflatex
pdf:	$(DOCDIR)/ug.pdf

$(DOCDIR)/ug.pdf: $(SRCFILES) $(PDFFIGS)
	$(ECHO) "Creating NAMD $(NAMD_VERSION) User's Guide ..." ; \
	$(CD) $(DOCDIR); \
	pdflatex $(MAIN_UG_FILE); \
	bibtex $(MAIN_UG_FILE); \
	pdflatex $(MAIN_UG_FILE); \
	pdflatex $(MAIN_UG_FILE); \

DISTILL = distill -compatlevel 2.1

# need to specify page size because distill can't figure it out (yet?)

figures/fmaOn.pdf: figures/fmaOn.eps
	$(DISTILL) -pagesize 301 174 pts < $^ > $@

figures/pairlistdist.pdf: figures/pairlistdist.eps
	$(DISTILL) -pagesize 240 240 pts < $^ > $@

figures/shifting.pdf: figures/shifting.eps
	$(DISTILL) -pagesize 301 165 pts < $^ > $@

figures/switching.pdf: figures/switching.eps
	$(DISTILL) -pagesize 310 165 pts < $^ > $@

figures/dual_top.pdf: figures/dual_top.eps
	$(DISTILL) -pagesize 584 220 pts < $^ > $@

$(DOCDIR)/ug.ps: $(SRCFILES)
	$(ECHO) "Creating NAMD $(NAMD_VERSION) User's Guide ..." ; \
	$(CD) $(DOCDIR); \
	latex $(MAIN_UG_FILE); \
	bibtex $(MAIN_UG_FILE); \
	latex $(MAIN_UG_FILE); \
	latex $(MAIN_UG_FILE); \
	dvips -t letter -o ug.ps ug.dvi

namd_version.tex:
	$(DELETE) namd_version.tex ; \
	echo '\\newcommand{\\NAMDVER}' "{$(NAMD_VERSION)}" > namd_version.tex

psfgen_macros.tex: ../psfgen/doc/ug_macros.tex
	$(DELETE) psfgen_macros.tex ; \
	grep PSFGEN ../psfgen/doc/ug_macros.tex > psfgen_macros.tex

ug_psfgen.tex: ../psfgen/doc/ug_psfgen.tex
	$(DELETE) ug_psfgen.tex ; \
	sed -e "20,$$ s/\\subsection{/\\subsubsection{/" \
	    -e "20,$$ s/\\section{/\\subsection{/" \
	../psfgen/doc/ug_psfgen.tex > ug_psfgen.tex

clean:
	$(DELETE) namd_version.tex ; \
	$(DELETE) psfgen_macros.tex ; \
	$(DELETE) ug_psfgen.tex ; \
	$(DELETE) $(DOCDIR)/ug.ps ; \
	$(DELETE) $(DOCDIR)/ug.aux ; \
	$(DELETE) $(DOCDIR)/ug.bbl ; \
	$(DELETE) $(DOCDIR)/ug.blg ; \
	$(DELETE) $(DOCDIR)/ug.dvi ; \
	$(DELETE) $(DOCDIR)/ug.lof ; \
	$(DELETE) $(DOCDIR)/ug.log ; \
	$(DELETE) $(DOCDIR)/ug.toc ; \
	$(DELETE) $(DOCDIR)/ug.pdf ; \
	$(DELETE) $(PDFFIGS)
