SYSTEMS = $(wildcard *x*x*_A*B*_bv?.??)
MODELS = $(foreach s,$(SYSTEMS),$(wildcard $(s)/temp*_exp*_den*_pop*))
SIMULATIONS_OUT = $(foreach m,$(MODELS),$(wildcard $(m)/*.out))
SIMULATIONS_CTF = $(foreach m,$(MODELS),$(wildcard $(m)/*.ctf))
SIMULATIONS_ENERGIES = $(foreach m,$(MODELS),$(wildcard $(m)/*.energy.npy.bz2))
SIMULATIONS_JPG = $(foreach m,$(MODELS),$(wildcard $(m)/*.jpg))

SYNC_REMOTE = ~/mnt/poly/scratch/
SYNC_OPTIONS = --verbose --progress --stats --human-readable --archive --compress --update
SYNC_EXCLUDES = --exclude *.ctf --exclude *.cga --exclude *.csa --exclude *.data-raw.npz
SYNC_DATA = /home/kml/data/

XVFBOPTS = "-screen 0 1280x1024x24"
PYCULGI = "xvfb-run -s $(XVFBOPTS) python-culgi"

.PHONY: default all
default: energy gallery

.PHONY: energy energy-total energy-field energy-coupling
energy: energy-total energy-field energy-coupling
energy-total: $(subst .npy.bz2,-total.png,$(SIMULATIONS_ENERGIES))
energy-field: $(subst .npy.bz2,-field.png,$(SIMULATIONS_ENERGIES))
energy-coupling: $(subst .npy.bz2,-coupling.png,$(SIMULATIONS_ENERGIES))

%.energy-total.png: %.energy.npy.bz2
	python plot-energy.py $< total save

%.energy-field.png: %.energy.npy.bz2
	python plot-energy.py $< field save

%.energy-coupling.png: %.energy.npy.bz2
	python plot-energy.py $< coupling save

.PHONY: convert
convert: $(subst .out,.data-raw.npz,$(SIMULATIONS_OUT))
%.data-raw.npz: %.csa %.cga %_Inst.ctf
	-"$(PYCULGI)" convert.py $(subst .csa,.out,$<)

.PHONY: analyze
analyze:  $(subst .out,.data-analyzed.npz,$(SIMULATIONS_OUT))
%.data-analyzed.npz: %.data-raw.npz
	python-culgi analyze.py $<

.PHONY: extract extract-energy extract-histograms
extract: extract-energy extract-histograms
extract-energy: $(subst .out,.energy.npy.bz2,$(SIMULATIONS_OUT))
extract-histograms: $(subst .out,.histograms.npy.bz2,$(SIMULATIONS_OUT))

%.energy.npy.bz2: %_Inst.ctf
	rm -rvf  $@
	-"$(PYCULGI)" extract.py $< && bzip2 $(subst .bz2,,$@)

%.histograms.npy.bz2: %.csa %.cga
	rm -rvf $@
	-"$(PYCULGI)" extract.py $(subst .csa,.out,$<) histograms && bzip2 $(subst .bz2,,$@)

.PHONY: copy
copy:
	rsync $(SYNC_OPTIONS) $(SYNC_EXCLUDES) $(SYNC_REMOTE)$(subst $(SYNC_DATA),,$(CURDIR))/*x*x*_A*B* .
	chmod 755 *x*x*_A*B*
	chmod 755 *x*x*_A*B*/*
	chmod 644 *x*x*_A*B*/*/*

.PHONY: avi
avi: $(subst .out,.avi,$(SIMULATIONS_OUT))

%.avi: %.csa %.cga
	-"$(PYCULGI)" replay.py $@ save xvfb && mencoder "mf://*.jpg" -mf fps=10 -o $@ -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600 && cp `ls *.jpg | tail -1` $(subst .avi,.jpg,$@)
	rm -rvf *.jpg

.PHONY: gallery
gallery: gallery.html
gallery.html: gallery.py $(SIMULATIONS_OUT) $(SIMULATIONS_JPG) */*/*.png
	python-culgi gallery.py > gallery.html
