all: lib zlib tinyxml
.PHONY: zlib tinyxml

lib:
	mkdir lib

zlib:
	make -C zlib
	cp zlib/build/libz.a lib/

clean_zlib:
	make -C zlib clean

tinyxml:
	make -C tinyxml
	cp tinyxml/build/libtinyxml.a lib/

clean_tinyxml:
	make -C tinyxml clean

clean: clean_zlib clean_tinyxml
	rm -f lib/libz.a lib/libtinyxml.a