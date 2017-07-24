#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Build supernova
#

VERSION=$(shell git describe --tags --always --dirty)


LIBPY=lib/python

LIBGO=$(shell pwd)/lib/go
LIBGO_BINS=reads/bucket_fastq_by_bc reads/sort_fastq_by_bc
export GOPATH=$(shell pwd)/lib/go:$(shell pwd)/tenkit/lib/go

LIBASSEMBLY=$(shell pwd)/lib/assembly
LIBASSEMBLY_BINS=CP DF MC TR ParseBarcodedFastqs FastFastbCount MakeFasta

LIBJEMALLOC=$(shell pwd)/lib/jemalloc
LIBJEMALLOC_BUILD=$(LIBJEMALLOC)/build
LIBJEMALLOC_CONF=$(LIBJEMALLOC_BUILD)/configure
LIBJEMALLOC_VERSION=3.6.0
LIBJEMALLOC_VERSIONED := $(LIBJEMALLOC)/$(LIBJEMALLOC_VERSION)

LIBTADA=$(shell pwd)/lib/tada
LIBTADA_BINS=tada
CARGO_DIR=$(shell pwd)


export DEPEND_FLAGS := -static

export CS=

.PHONY: $(LIBTADA_BINS) $(LIBASSEMBLY_BINS) $(LIBGO_BINS) test clean jemalloc

#
# Targets for development builds.
#
all: $(LIBTADA_BINS) $(LIBGO_BINS) $(LIBJEMALLOC_VERSIONED) $(LIBASSEMBLY_BINS) tenkit-all

tenkit-all:
	make -C tenkit all

$(LIBJEMALLOC_VERSIONED): $(LIBJEMALLOC_CONF)
	make -C $(LIBJEMALLOC_BUILD)
	make -C $(LIBJEMALLOC_BUILD) install_lib install_bin install_include DESTDIR=..

$(LIBJEMALLOC_CONF):
	cd $(LIBJEMALLOC_BUILD); \
	./autogen.sh ;\
	CPPFLAGS="-Wno-unused -Wno-maybe-uninitialized" ./configure --prefix=/$(LIBJEMALLOC_VERSION) --with-jemalloc-prefix=je_

jemalloc: $(LIBJEMALLOC_VERSIONED)

lib/bin:
	mkdir lib/bin

$(LIBGO_BINS): lib/bin
	go install -ldflags "-X $@.__VERSION__ '$(VERSION)'" $@
	cp lib/go/bin/* lib/bin

$(LIBASSEMBLY_BINS): lib/bin $(LIBASSEMBLY) jemalloc
	$(MAKE) -j16 -C $(LIBASSEMBLY) STATIC=yes $@
	cp $(LIBASSEMBLY)/bin/$@ lib/bin

$(LIBTADA_BINS): lib/bin $(LIBTADA)
	cd $(LIBTADA); CC=gcc CXX=g++ CARGO_HOME=$(CARGO_DIR)/.cargo cargo build --release
	cp $(LIBTADA)/target/release/$@ lib/bin


clean: 
	rm -rf $(LIBGO)/bin
	rm -rf $(LIBGO)/pkg
	rm -rf $(LIBASSEMBLY)/bin
	rm -rf $(CARGO_DIR)/.cargo
	rm -rf $(TADAPATH)/target
	rm -rf lib/bin
	rm -rf $(LIBJEMALLOC_VERSIONED)
	rm -rf $(LIBJEMALLOC_CONF)
	make -C $(LIBASSEMBLY) cleanAll


