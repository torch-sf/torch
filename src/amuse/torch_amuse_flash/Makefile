ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
FLASH_VERSION := FLASH4.6.2
export FLASH_VERSION

DEPFLAGS += $(FORSOCKETS_CFLAGS) $(OPENMP_FFLAGS) $(HDF5_FFLAGS)

FCFLAGS += -fallow-argument-mismatch $(DEPFLAGS)

LDFLAGS += $(HDF5_LDFLAGS)

LDLIBS += -lm $(FORSOCKETS_LIBS) $(HDF5_FLIBS)


CODELIB = src/libflash.a


.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all AR=$(AR) RANLIB=$(RANLIB)


# Building the workers
flash_worker.f90: interface.py
	amusifier --type=f90 interface.py FlashInterface -o $@

flash_worker.o: flash_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

flash_worker: base_grid_interface.o interface.o flash_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


%.o: %.F90 $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Isrc/$(FLASH_VERSION)/object $<


# Which packages contain which workers?
torch-amuse-flash_contains: flash_worker


# Building and installing packages
develop-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install -e packages/$*

install-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install packages/$*

package-%: %_contains
	python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} packages/$*

test-%:
	cd packages/$* && pytest


# Cleaning up
.PHONY: clean
clean:
	$(MAKE) -C src clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
