DUMMY := $(shell $(MAKE) DEBUGGING=$(DEBUGGING) PROFILE=$(PROFILE) GENERIC=$(GENERIC) -fMakefile_init)

.PHONY: all
all: install_scripts

include Makefile_scripts
include Makefile_rules
include Makefile_deps
include Makefile_special
$(LIBNAM): $(LIBOBJECTS)
