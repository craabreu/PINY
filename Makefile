SYSTEMS=$(filter-out make_defs, $(notdir $(wildcard ./compile/*)))
CLEANERS=$(addsuffix -clean,${SYSTEMS})
DEFAULT=linux

.PHONY: all $(SYSTEMS) $(CLEANERS) clean-all

all: $(DEFAULT)

clean: $(DEFAULT)-clean

clean-all: $(CLEANERS)

$(SYSTEMS):
	make -C ./compile/$@

$(CLEANERS):
	make -C ./compile/$(@:%-clean=%) clean
