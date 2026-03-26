CODELIB = libflash.a

all: $(CODELIB) 

$(CODELIB): $(CODEOBJS)
	make -C $(FLASH_VERSION)/object
	rm -f $@
	$(AR) -r $@ $(FLASH_VERSION)/object/*.o
	$(AR) -d $@ Flash.o UNIT_TEST.o
	$(RANLIB) $@

clean_codelib:
	rm -f $(CODELIB)

clean: clean_codelib
	make -C $(FLASH_VERSION)/object clean

distclean: clean_codelib
	make -C $(FLASH_VERSION)/object distclean

