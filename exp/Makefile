# SEM images and outputs.
SEM_JPGS = $(wildcard sem/*/C*/*.jpg)
SEM_RDFS = $(subst sem/,sem-analyzed/,$(subst .jpg,-rdf.png,$(SEM_JPGS)))

# Frame buffer parameters (xvfb).
XVFBOPTS = "-screen 11 1280x1024x24"
PYTHON = "xvfb-run -s $(XVFBOPTS) python-culgi"

.PHONY: default all benchmark sem sem-comp
default: all
all: benchmark sem sem-comp

benchmark: benchmark/F20-atom_res.preview-rdf.png
benchmark/%-rdf.png: benchmark/%.png
	"$(PYTHON)" benchmark.py $< save && bzip2 benchmark/$**.npy

sem: $(SEM_RDFS)
sem-analyzed/%-rdf.png: sem/%.jpg
	"$(PYTHON)" sem.py $< save && bzip2 $(subst -rdf.png,-*.npy,$@)

#sem-comp: $(foreach c,$(SEM_CONC),sem-rdf.$(c).png)

#sem/rdf/%.png: sem/%.tif
#	python analyze.py $< save

#sem-rdf.%.png: sem/%.*.tif
#	python analyze.py comp $^
