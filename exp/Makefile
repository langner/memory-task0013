# Frame buffer parameters (xvfb)
XVFBOPTS = "-screen 11 1280x1024x24"
PYTHON = "xvfb-run -s $(XVFBOPTS) python-culgi"

.PHONY: default all benchmark sem sem-comp
default: all
all: benchmark sem sem-comp

# A becnhmark
BENCHMARK = benchmark/F20-atom_res.preview
benchmark: $(BENCHMARK)-rdf.png
$(BENCHMARK)-rdf.png: $(BENCHMARK).png benchmark.py common.py
	rm -rvf $@
	"$(PYTHON)" benchmark.py $< save
	-bzip2 benchmark/$**.npy
	rm -rvf benchmark/$**.npy
	sleep 2

# SEM images
SEM_JPGS = $(wildcard sem/*/C*/*.jpg)
SEM_NPS = $(subst sem/,sem-analyzed/,$(subst .jpg,-nps.png,$(SEM_JPGS)))
SEM_RDF = $(subst sem/,sem-analyzed/,$(subst .jpg,-rdf.png,$(SEM_JPGS)))
sem: $(SEM_NPS) $(SEM_RDF)
sem-analyzed/%-nps.png sem-analyzed/%-rdf.png: sem/%.jpg
	"$(PYTHON)" sem.py $< save
	-bzip2 $**.npy
	rm -rvf $**.npy
	sleep 2

.PHONY: clean
clean:
	rm -rvf benchmark/*.npy* benchmark/*-{coms,filtered,nps,rdf,threshold}.png
	rm -rvf sem-analyzed/*/*/*.{npy,bz2,png}
