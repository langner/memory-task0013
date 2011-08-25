SYSTEMS = $(wildcard *x*x*_A*B*_bv?.??)
MODELS = $(foreach s,$(SYSTEMS),$(wildcard $(s)/temp*_exp*_den*_pop*))
SIMULATIONS_OUT = $(foreach m,$(MODELS),$(wildcard $(m)/*.out))
SIMULATIONS_CTF = $(foreach m,$(MODELS),$(wildcard $(m)/*.ctf))
SIMULATIONS_ENERGIES = $(foreach m,$(MODELS),$(wildcard $(m)/*.energy.npy.bz2))

SYNC_REMOTE = ~/mnt/poly/scratch/
SYNC_OPTIONS = --verbose --progress --stats --human-readable --archive --compress --update
SYNC_EXCLUDES = --exclude *.ctf --exclude *.cga --exclude *.csa
SYNC_DATA = /home/kml/data/

XVFBOPTS = "-screen 0 1280x1024x24"
PYCULGI = "xvfb-run -s $(XVFBOPTS) python-culgi"

.PHONY: default all
default: energy gallery

.PHONY: energy
energy: $(subst .npy.bz2,.png,$(SIMULATIONS_ENERGIES))

%.energy.png: %.energy.npy.bz2 plot-energy.py
	python plot-energy.py $< save

.PHONY: extract extract-energy extract-csa extract-cga
extract: extract-energy # extract-csa extract-cga
extract-energy: $(subst .out,.energy.npy.bz2,$(SIMULATIONS_OUT))
extract-csa: $(subst .out,.csa.npy.bz2,$(SIMULATIONS_OUT))
extract-cga: $(subst .out,.cga.npy.bz2,$(SIMULATIONS_OUT))

%.energy.npy.bz2: %_Inst.ctf
	rm -rvf  $@
	-"$(PYCULGI)" extract.py $< && bzip2 $(subst .bz2,,$@)

%.csa.npy.bz2: %.csa
	rm -rvf $@
	-python-culgi extract.py archives $< && bzip2 $(subst .bz2,,$@)

%.cga.npy.bz2: %.cga
	rm -rvf $@
	-python-culgi extract.py archives $< && bzip2 $(subst .bz2,,$@)

#%.histograms.npy.bz2: %.csa %.cga
#	rm -rvf  $@
#	-python-culgi extract.py $(subst .csa,.out,$<) archive && bzip2 $(subst .bz2,,$@)

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
gallery.html: gallery.py $(SIMULATIONS_OUT)
	python-culgi gallery.py > gallery.html
