# Frame buffer parameters (xvfb)
XVFBOPTS = "-screen 11 1280x1024x24"
PYTHON = "xvfb-run -s $(XVFBOPTS) python-culgi"

.PHONY: default all benchmark sem sem-comp
default: all
all: benchmark sem sem-comp article0015

# A becnhmark.
BENCHMARKS = benchmark/F20-atom_res.preview benchmark/a06fig02a benchmark/r6 benchmark/r7
BENCHMARKS_RDF = $(foreach b,$(BENCHMARKS),$(b)-rdf.png)
benchmark: $(BENCHMARKS_RDF)
%-rdf.png: %.png benchmark.py common.py
	rm -rvf $@
	"$(PYTHON)" benchmark.py $< save
	-bzip2 benchmark/$**.npy
	rm -rvf benchmark/$**.npy
	sleep 2

# SEM images.
SEM_JPGS = $(wildcard sem/*/C*/*.jpg)
SEM_NPS = $(subst sem/,sem-analyzed/,$(subst .jpg,-nps.png,$(SEM_JPGS)))
SEM_RDF = $(subst sem/,sem-analyzed/,$(subst .jpg,-rdf.png,$(SEM_JPGS)))
sem: $(SEM_NPS) $(SEM_RDF)
sem-analyzed/%-nps.png sem-analyzed/%-rdf.png: sem/%.jpg
	"$(PYTHON)" sem.py $< save
	rm -rvf $**.npy.bz2
	-bzip2 sem-analyzed/$**.npy
	rm -rvf $**.npy
	sleep 2

# Plots for article0015.
.PHONY: article0015
article0015: article0015-fig2.png article0015-fig3.png article0015-fig4-single.png article0015-fig4.png
article0015-fig4-single: article0015.py
	python article0015.py fig4 single
article0015-%.png: article0015.py
	python article0015.py $*

.PHONY: clean
clean:
	rm -rvf benchmark/*.npy* benchmark/*-{coms,filtered,nps,rdf,threshold}.png
	rm -rvf sem-analyzed/*/*/*.{npy,bz2,png}
