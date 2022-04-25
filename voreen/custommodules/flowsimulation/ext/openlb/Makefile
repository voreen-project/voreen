###########################################################################
## installing the configuration makefile

ifdef IN_NIX_SHELL
  ifndef OPENLB_CONFIG
    OPENLB_CONFIG = config.mk
  endif
  MAKEFILE_INC_TEMP := $(shell cp -n ${OPENLB_CONFIG} config.mk)
endif

CLEAN_TARGETS :=

###########################################################################
## Embedded dependencies (optional)

all: dependencies

dependencies:
	$(MAKE) -C external

clean-dependencies:
	make -C external clean

CLEAN_TARGETS += clean-dependencies

###########################################################################
## code generation (CSE)

CSE_GENERATORS := $(shell find src/ -type f -name '*.cse.h.template')
CSE_GENERATEES := $(patsubst %.cse.h.template, %.cse.h, $(CSE_GENERATORS))

%.cse.h: %.cse.h.template
	python codegen/cse.py $< $@

cse: $(CSE_GENERATEES)

###########################################################################
## Examples

EXAMPLES := $(dir $(shell find examples -name 'Makefile'))

$(EXAMPLES):
	$(MAKE) -C $@ onlysample

.PHONY: $(EXAMPLES)

samples: dependencies $(EXAMPLES)

###########################################################################
## Cleaning

clean: $(CLEAN_TARGETS)

###########################################################################
## User guide documentation

userguide:
	@cd doc/userGuide/; \
	latexmk -pdf -silent -f olb-ug.tex

###########################################################################
## Doxygen documentation

doxygen:
	doxygen doc/DoxygenConfig