# Normally do plots and gallery, perhaps copy.
.PHONY: default all
default: plot gallery report exp
all: copy default

include $(wildcard ../common.mk)

# ############
# Output files
# ############

# Simulation output files
PHASES = $(wildcard phase?)
SYSTEMS = $(foreach p,$(PHASES),$(wildcard $(p)/*x*x*_A*B*_bv?.??))
MODELS = $(foreach s,$(SYSTEMS),$(wildcard $(s)/temp*_exp*_den*_pop*))
SIMS = $(foreach m,$(MODELS),$(wildcard $(m)/k*_nchi*))
SIMS_OUT = $(foreach s,$(SIMS),$(wildcard $(s)/*.out))
SIMS_OUT_NEW = $(subst .cga,.out,$(foreach s,$(SIMS),$(wildcard $(s)/*.cga)))
SIMS_JPG = $(foreach m,$(SIMS),$(wildcard $(m)/*.jpg))
SIMS_AVI = $(foreach m,$(SIMS),$(wildcard $(m)/*.avi))
SIMS_ENERGY = $(foreach m,$(SIMS),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_HIST_FIELD = $(foreach m,$(SIMS),$(wildcard $(m)/*.hist-field.npy.bz2))
SIMS_HIST_RADIAL = $(foreach m,$(SIMS),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RESIDUAL = $(foreach m,$(SIMS),$(wildcard $(m)/*.hist-residual.npy.bz2))

# Plots to be generated.
PLOTS_ENERGY_TOTAL = $(subst .energy.npy.bz2,.energy-total.png,$(SIMS_ENERGY))
PLOTS_ENERGY_FIELD = $(subst .energy.npy.bz2,.energy-field.png,$(SIMS_ENERGY))
PLOTS_HIST_FIELD_TOTAL = $(subst .hist-field.npy.bz2,.hist-field-total.png,$(SIMS_HIST_FIELD))
PLOTS_HIST_FIELD_ORDER = $(subst .hist-field.npy.bz2,.hist-field-order.png,$(SIMS_HIST_FIELD))
PLOTS_HIST_RADIAL = $(subst .hist-radial.npy.bz2,.hist-radial.png,$(SIMS_HIST_RADIAL))
PLOTS_HIST_RADIAL_ZOOM = $(subst .hist-radial.npy.bz2,.hist-radial-zoom.png,$(SIMS_HIST_RADIAL))
PLOTS_HIST_RESIDUAL_TOTAL = $(subst .hist-residual.npy.bz2,.hist-residual-total.png,$(SIMS_HIST_RESIDUAL))
PLOTS_HIST_RESIDUAL_ORDER = $(subst .hist-residual.npy.bz2,.hist-residual-order.png,$(SIMS_HIST_RESIDUAL))
PLOTS_ENERGY = $(PLOTS_ENERGY_TOTAL) $(PLOTS_ENERGY_FIELD)
PLOTS_HIST_FIELD = $(PLOTS_HIST_FIELD_TOTAL) $(PLOTS_HIST_FIELD_ORDER)
PLOTS_HIST_RESIDUAL = $(PLOTS_HIST_RESIDUAL_TOTAL) $(PLOTS_HIST_RESIDUAL_ORDER)

# Synchronization parameters (rsync)
SYNC_REMOTE = poly:scratch/
SYNC_OPTIONS = --verbose --progress --stats --human-readable --archive --compress --update
SYNC_EXCLUDES = --exclude *.ctf* --exclude *.cga* --exclude *.csa*
SYNC_DATA = /home/kml/data/

# Frame buffer parameters (xvfb)
XVFBOPTS = "-screen 0 1280x1024x24"
PYCULGI = "xvfb-run -s $(XVFBOPTS) python-culgi"
PYTHON = "xvfb-run -n $(shell echo $$RANDOM) -s $(XVFBOPTS) python"

# #############
# Local targets
# #############

# Data is analyzed on cluster, so synchronize local files with that
.PHONY: copy
copy:
	rsync $(SYNC_OPTIONS) $(SYNC_EXCLUDES) $(SYNC_REMOTE)$(subst $(SYNC_DATA),,$(CURDIR))/phase? .
	chmod 755 phase?
	chmod 755 phase?/*x*x*_A*B*
	chmod 755 phase?/*x*x*_A*B*/*
	chmod 755 phase?/*x*x*_A*B*/*/*
	chmod 644 phase?/*x*x*_A*B*/*/*/*

# Generate the plots based on analyzed data
# Don't depend on coupling and offset plots, since neat systems don't get those
.PHONY: plot plot-energy plot-offsets plot-hist-field plot-hist-radial plot-hist-residual
plot: plot-energy plot-offsets plot-hist-field plot-hist-radial plot-hist-residual
plot-energy: $(PLOTS_ENERGY)
%.energy-total.png %.energy-field.png: %.energy.npy.bz2
	"$(PYTHON)" plot.py $< total save
	"$(PYTHON)" plot.py $< field save
	"$(PYTHON)" plot.py $< coupl save
	"$(PYTHON)" plot.py $< offsets save
plot-hist-field: $(PLOTS_HIST_FIELD)
%.hist-field-total.png %.hist-field-order.png: %.hist-field.npy.bz2
	"$(PYTHON)" plot.py $< total save
	"$(PYTHON)" plot.py $< order save
plot-hist-radial: $(PLOTS_HIST_RADIAL) $(PLOTS_HIST_RADIAL_ZOOM)
%.hist-radial.png %.hist-radial-zoom.png: %.hist-radial.npy.bz2
	"$(PYTHON)" plot.py $< save
	"$(PYTHON)" plot.py $< zoom save
plot-hist-residual: $(PLOTS_HIST_RESIDUAL)
%.hist-residual-total.png %.hist-residual-order.png: %.hist-residual.npy.bz2
	"$(PYTHON)" plot.py $< total save
	"$(PYTHON)" plot.py $< order save

# Generate the galleries if any key output files changed
.PHONY: gallery
gallery: gallery.html $(foreach p,$(PHASES),$(p)/gallery.html) exp/gallery.html
gallery.html: gallery.py
	python-culgi gallery.py main > gallery.html && rm -rvf Culgi.log
phase%/gallery.html: phase%/*/*/*/*.out phase%/*/*/*/*.jpg phase%/*/*/*/*.avi phase%/*/*/*/*.png gallery.py
	python-culgi gallery.py $* > phase$*/gallery.html && rm -rvf Culgi.log
exp/gallery.html: gallery.py exp/sem-analyzed/*/*/*.png
	python-culgi gallery.py exp > $@ && rm -rvf Culgi.log

# Generate report.
.PHONY: report
report: task.pdf

# Analyze experimental data.
.PHONY: exp
exp:
	$(MAKE) -c exp

# #############################
# Remote targets (run manually)
# #############################

# A safe sequence of operations on the remote target would be:
#  1. avi
#  2. convert
#  3. analyze
# The snapshot/movie target comes first, because after convert
#  the archives will be compressed to save space, which is not
#  a problem for analysis, but for Culgi yes
# Note that currently I need to do avi on a different machine
#  via sshfs, because no encoder is installed on poly

# Generate movies and snapshots for simulations based on Culgi outputs
.PHONY: avi
avi: $(subst .out,.avi,$(SIMS_OUT_NEW))
%.avi: %.csa %.cga
	-"$(PYCULGI)" replay.py $@ save xvfb && mencoder "mf://*.jpg" -mf fps=10 -o $@ -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600 && cp `ls *.jpg | head -1` $(subst .avi,.first.jpg,$@) && cp `ls *.jpg | tail -1` $(subst .avi,.jpg,$@)
	rm -rvf *.jpg

# Convert Culgi output (cga/csa/ctf) to NumPy compressed archives (npz)
# Also compress both targets and dependencies after the conversion
.PHONY: convert
convert: $(subst .out,.ctf.npy.gz,$(SIMS_OUT_NEW)) $(subst .out,.csa.npy.gz,$(SIMS_OUT_NEW)) $(subst .out,.cga.npy.gz,$(SIMS_OUT_NEW))
%.ctf.npy.gz %.csa.npy.gz %.cga.npy.gz: %.out %_Inst.ctf %.csa %.cga
	-"$(PYCULGI)" convert.py $< && gzip $*{_Inst.ctf,.csa,.cga} && gzip $*.{ctf,csa,cga}.npy

# Analyze the data (npy format, compressed)
# The resulting archives should be much smaller than the raw data
# Do not depend explicitely on csa files, which are missing for neat systems (energy and fields should suffice)
.PHONY: analyze
analyze: $(subst .out,.energy.npy.bz2,$(SIMS_OUT)) $(subst .out,.hist-field.npy.bz2,$(SIMS_OUT))
%.energy.npy.bz2 %.hist-field.npy.bz2: %.out %.ctf.npy.gz %.cga.npy.gz
	-python-culgi analyze.py $< energy histograms && bzip2 $*.energy*.npy $*.hist*.npy
