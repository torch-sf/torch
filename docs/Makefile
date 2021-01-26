SUBDIRS = git_workflow sinks torch_manual quickstart

all:
	cd git_workflow && $(MAKE)
	cd sinks && $(MAKE)
	cd torch_manual && $(MAKE)
	cd quickstart && $(MAKE)


.PHONY: clean

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done

