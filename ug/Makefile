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

### latex the documentation
doc:	$(DOCDIR)/ug.ps 

$(DOCDIR)/ug.ps: 
	$(ECHO) "Creating User's Guide ..." ; \
	$(CD) $(DOCDIR); \
	latex $(MAIN_UG_FILE); \
	bibtex $(MAIN_UG_FILE); \
	latex $(MAIN_UG_FILE); \
	latex $(MAIN_UG_FILE); \
	dvips -o ug.ps ug.dvi


clean:
	$(DELETE) $(DOCDIR)/ug.ps ; \
	$(DELETE) $(DOCDIR)/ug.aux ; \
	$(DELETE) $(DOCDIR)/ug.bbl ; \
	$(DELETE) $(DOCDIR)/ug.blg ; \
	$(DELETE) $(DOCDIR)/ug.lof ; \
	$(DELETE) $(DOCDIR)/ug.log ; \
	$(DELETE) $(DOCDIR)/ug.toc

