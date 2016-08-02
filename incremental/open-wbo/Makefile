VERSION    = core
SOLVERNAME = "laser"
SOLVERDIR  = laser
NSPACE     = Minisat
# THE REMAINING OF THE MAKEFILE SHOULD BE LEFT UNCHANGED
EXEC       = open-wbo
DEPDIR     = mtl utils core 
DEPDIR     +=  ../../encodings ../../algorithms
MROOT      = $(PWD)/solvers/$(SOLVERDIR)
CFLAGS     = -Wall -Wno-parentheses -DNSPACE=$(NSPACE) -DSOLVERNAME=$(SOLVERNAME) -DVERSION=$(VERSION)
ifeq ($(VERSION),simp)
DEPDIR     += simp
CFLAGS     += -DSIMP=1 
endif
include $(MROOT)/mtl/template.mk
