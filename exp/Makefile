SHELL=/bin/bash

# Frame buffer parameters (xvfb)
XVFBOPTS = "-screen 0 1280x1024x24"
PYCULGI = python-culgi-6.1.0-kml
PYCULGI_FB = "xvfb-run -n $(shell echo $$RANDOM) -s $(XVFBOPTS) $(PYCULGI)"

.PHONY: default all benchmark sem homo
default: all
all: benchmark sem analysis

# A few benchmarks
BENCHMARKS = benchmark/F20-atom_res.preview benchmark/a06fig02a benchmark/r6 benchmark/r7
BENCHMARKS_RDF = $(foreach b,$(BENCHMARKS),$(b)-rdf.png)
benchmark: $(BENCHMARKS_RDF)
%-rdf.png: %.png benchmark.py common.py
	rm -rvf $@
	"$(PYCULGI_FB)" benchmark.py $< save
	-bzip2 benchmark/$**.npy
	rm -rvf benchmark/$**.npy
	sleep 2

# SEM images
SEM_JPGS = $(wildcard sem/*/C*/*.jpg) 
SEM_RDF = $(subst sem/,sem-analyzed/,$(subst .jpg,-rdf.png,$(SEM_JPGS)))
sem: $(SEM_RDF)
sem-analyzed/%-nps.png sem-analyzed/%-rdf.png: sem/%.jpg
	"$(PYCULGI_FB)" sem.py $< save
	rm -rvf $**.npy.bz2
	-bzip2 sem-analyzed/$**.npy
	rm -rvf $**.npy
	sleep 2

# Homopolymer images
HOMO_JPGS = $(wildcard homo/*/C*/*.jpg) 
HOMO_RDF = $(subst homo/,homo-analyzed/,$(subst .jpg,-rdf.png,$(HOMO_JPGS)))
homo: $(HOMO_RDF)
homo-analyzed/%-nps.png homo-analyzed/%-rdf.png: homo/%.jpg
	"$(PYCULGI_FB)" sem.py $< save
	rm -rvf $**.npy.bz2
	-bzip2 homo-analyzed/$**.npy
	rm -rvf $**.npy
	sleep 2

# The a .... finish this

.PHONY: clean
clean:
	rm -rvf benchmark/*.npy* benchmark/*-{coms,filtered,nps,rdf,threshold}.png
	rm -rvf sem-analyzed/*/*/*.{npy,bz2,png}