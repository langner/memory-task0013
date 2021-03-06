SHELL=/bin/bash

# By default do plots and gallery, perhaps copy.
.PHONY: default all
default: plot gallery report exp
all: copy default

include $(wildcard ../common.mk)

AVIDIR := $(shell mktemp -d -u)

# ############
# Output files
# ############

# Simulation output files
PHASES = $(wildcard phase?) $(wildcard phase??)
SYSTEMS = $(foreach p,$(PHASES),$(wildcard $(p)/*x*x*_A*B*_bv*))
SYSTEMS8 = $(wildcard phase8/*x*x*_A*B*_bv?.??)
SYSTEMS9 = $(wildcard phase9/*x*x*_A*B*_bv?.??)
SYSTEMS10 = $(wildcard phase1?/*x*x*_A*B*_bv*)
SYSTEMS17 = $(wildcard phase1[7-9]/*x*x*_A*B*_bv*)
SYSTEMS20 = $(wildcard phase2?/*x*x*_A*B*_bv*)
MODELS = $(foreach s,$(SYSTEMS),$(wildcard $(s)/temp*_exp*_den*_pop*))
MODELS8 = $(foreach s,$(SYSTEMS8),$(wildcard $(s)/temp*_exp*_den*_pop*))
MODELS9 = $(foreach s,$(SYSTEMS9),$(wildcard $(s)/temp*_exp*_den*_pop*))
MODELS10 = $(foreach s,$(SYSTEMS10),$(wildcard $(s)/temp*_exp*_den*_pop*))
MODELS17 = $(foreach s,$(SYSTEMS17),$(wildcard $(s)/temp*_exp*_den*_pop*))
MODELS20 = $(foreach s,$(SYSTEMS20),$(wildcard $(s)/temp*_exp*_den*_pop*))
SIMS = $(foreach m,$(MODELS),$(wildcard $(m)/k*_nchi*))
SIMS8 = $(foreach m,$(MODELS8),$(wildcard $(m)/k*_nchi*))
SIMS9 = $(foreach m,$(MODELS9),$(wildcard $(m)/k*_nchi*))
SIMS10 = $(foreach m,$(MODELS10),$(wildcard $(m)/k*_nchi*))
SIMS17 = $(foreach m,$(MODELS17),$(wildcard $(m)/k*_nchi*))
SIMS20 = $(foreach m,$(MODELS20),$(wildcard $(m)/k*_nchi*))
SIMS_OUT = $(foreach s,$(SIMS),$(wildcard $(s)/*.out))

# More output files
SIMS_OUT_NEW = $(subst .cga,.out,$(foreach s,$(SIMS),$(wildcard $(s)/*.cga)))
SIMS_JPG = $(foreach m,$(SIMS),$(wildcard $(m)/*.jpg))
SIMS_AVI = $(foreach m,$(SIMS),$(wildcard $(m)/*.avi))
SIMS_ENERGY = $(foreach m,$(SIMS),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_ENERGY8 = $(foreach m,$(SIMS8),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_ENERGY9 = $(foreach m,$(SIMS9),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_ENERGY10 = $(foreach m,$(SIMS10),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_ENERGY17 = $(foreach m,$(SIMS17),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_ENERGY20 = $(foreach m,$(SIMS20),$(wildcard $(m)/*.energy.npy.bz2))
SIMS_HIST_FIELD = $(foreach m,$(SIMS),$(wildcard $(m)/*.hist-field.npy.bz2))
SIMS_HIST_RADIAL = $(foreach m,$(SIMS),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RADIAL8 = $(foreach m,$(SIMS8),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RADIAL9 = $(foreach m,$(SIMS9),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RADIAL10 = $(foreach m,$(SIMS10),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RADIAL17 = $(foreach m,$(SIMS17),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RADIAL20 = $(foreach m,$(SIMS20),$(wildcard $(m)/*.hist-radial.npy.bz2))
SIMS_HIST_RESIDUAL = $(foreach m,$(SIMS),$(wildcard $(m)/*.hist-residual.npy.bz2))
SIMS_HIST_RESIDUAL8 = $(foreach m,$(SIMS8),$(wildcard $(m)/*.hist-residual.npy.bz2))
SIMS_HIST_RESIDUAL9 = $(foreach m,$(SIMS9),$(wildcard $(m)/*.hist-residual.npy.bz2))
SIMS_HIST_RESIDUAL10 = $(foreach m,$(SIMS10),$(wildcard $(m)/*.hist-residual.npy.bz2))
SIMS_HIST_RESIDUAL17 = $(foreach m,$(SIMS17),$(wildcard $(m)/*.hist-residual.npy.bz2))
SIMS_HIST_RESIDUAL20 = $(foreach m,$(SIMS20),$(wildcard $(m)/*.hist-residual.npy.bz2))
SIMS_HIST_ANGLES = $(foreach m,$(SIMS9) $(SIMS10) $(SIMS20),$(wildcard $(m)/*.hist-ang.npy.bz2))

# Plots to be generated
PLOTS_ENERGY_TOTAL = $(subst .energy.npy.bz2,.energy-total.png,$(SIMS_ENERGY))
PLOTS_ENERGY_FIELD = $(subst .energy.npy.bz2,.energy-field.png,$(SIMS_ENERGY))
PLOTS_ENERGY_COUPL = $(subst .hist-residual.npy.bz2,.energy-coupl.png,$(SIMS_HIST_RESIDUAL))
PLOTS_OFFSETS_RAD = $(subst .hist-radial.npy.bz2,.offsets.png,$(SIMS_HIST_RADIAL8) $(SIMS_HIST_RADIAL9) $(SIMS_HIST_RADIAL10) $(SIMS_HIST_RADIAL20))
PLOTS_OFFSETS_ANG = $(subst .hist-radial.npy.bz2,.offsets-ang.png,$(SIMS_HIST_RADIAL9) $(SIMS_HIST_RADIAL10) $(SIMS_HIST_RADIAL20))
PLOTS_INTERFACE = $(subst .hist-residual.npy.bz2,.interface.png,$(SIMS_HIST_RESIDUAL17) $(SIMS_HIST_RESIDUAL20))
PLOTS_HIST_FIELD_TOTAL = $(subst .hist-field.npy.bz2,.hist-field-total.png,$(SIMS_HIST_FIELD))
PLOTS_HIST_FIELD_ORDER = $(subst .hist-field.npy.bz2,.hist-field-order.png,$(SIMS_HIST_FIELD))
PLOTS_HIST_RADIAL = $(subst .hist-radial.npy.bz2,.hist-radial.png,$(SIMS_HIST_RADIAL))
PLOTS_HIST_RADIAL_ZOOM = $(subst .hist-radial.npy.bz2,.hist-radial.zoom.png,$(SIMS_HIST_RADIAL))
PLOTS_HIST_RADIAL_SHELL = $(subst .hist-radial.npy.bz2,.hist-radial-shell.png,$(SIMS_HIST_RADIAL8) $(SIMS_HIST_RADIAL9) $(SIMS_HIST_RADIAL10) $(SIMS_HIST_RADIAL20))
PLOTS_HIST_RADIAL_SHELL_ZOOM = $(subst .hist-radial.npy.bz2,.hist-radial-shell.zoom.png,$(SIMS_HIST_RADIAL8) $(SIMS_HIST_RADIAL9) $(SIMS_HIST_RADIAL10) $(SIMS_HIST_RADIAL20))
PLOTS_HIST_RES_TOTAL = $(subst .hist-residual.npy.bz2,.hist-residual-total.png,$(SIMS_HIST_RESIDUAL))
PLOTS_HIST_RES_TOTAL_SHELL = $(subst .hist-residual.npy.bz2,.hist-residual-total-shell.png,$(SIMS_HIST_RESIDUAL8) $(SIMS_HIST_RESIDUAL9) $(SIMS_HIST_RESIDUAL10) $(SIMS_HIST_RESIDUAL20))
PLOTS_HIST_RES_ORDER = $(subst .hist-residual.npy.bz2,.hist-residual-order.png,$(SIMS_HIST_RESIDUAL))
PLOTS_HIST_RES_ORDER_SHELL = $(subst .hist-residual.npy.bz2,.hist-residual-order-shell.png,$(SIMS_HIST_RESIDUAL8) $(SIMS_HIST_RESIDUAL9) $(SIMS_HIST_RESIDUAL10) $(SIMS_HIST_RESIDUAL20))
PLOTS_ENERGY = $(PLOTS_ENERGY_TOTAL) $(PLOTS_ENERGY_FIELD) $(PLOTS_ENERGY_COUPL)
PLOTS_OFFSETS = $(PLOTS_OFFSETS_RAD) $(PLOTS_OFFSETS_ANG)
PLOTS_HIST_FIELD = $(PLOTS_HIST_FIELD_TOTAL) $(PLOTS_HIST_FIELD_ORDER)
PLOTS_HIST_RADIALS = $(PLOTS_HIST_RADIAL) $(PLOTS_HIST_RADIAL_ZOOM) $(PLOTS_HIST_RADIAL_SHELL) $(PLOTS_HIST_RADIAL_SHELL_ZOOM)
PLOTS_HIST_RES = $(PLOTS_HIST_RES_TOTAL) $(PLOTS_HIST_RES_TOTAL_SHELL) $(PLOTS_HIST_RES_ORDER) $(PLOTS_HIST_RES_ORDER_SHELL)
PLOTS_HIST_ANGLES = $(subst .hist-ang.npy.bz2,.hist-ang.png,$(SIMS_HIST_ANGLES))

# Synchronization parameters (for rsync)
SYNC_REMOTE = poly:scratch/
SYNC_OPTIONS = --verbose --progress --stats --human-readable --archive --compress --update
SYNC_EXCLUDES = --exclude *.ctf* --exclude *.cga* --exclude *.csa*
SYNC_DATA = /home/kml/data/

# Frame buffer parameters (for xvfb)
XVFBOPTS = "-screen 0 1280x1024x24"
PYCULGI = python-culgi-6.1.0-kml
PYCULGI_FB = "xvfb-run -n $(shell echo $$RANDOM) -s $(XVFBOPTS) $(PYCULGI)"
PYTHON = python
PYTHON_FB = "xvfb-run -n $(shell echo $$RANDOM) -s $(XVFBOPTS) $(PYTHON)"

# #############
# Local targets
# #############

# Synchronize local files with remote location.
.PHONY: copy
copy:
	rsync $(SYNC_OPTIONS) $(SYNC_EXCLUDES) $(SYNC_REMOTE)$(subst $(SYNC_DATA),,$(CURDIR))/phase* .
	chmod 755 phase*
	chmod 755 phase*/*x*x*_A*B*
	chmod 755 phase*/*x*x*_A*B*/*
	chmod 755 phase*/*x*x*_A*B*/*/*
	chmod 644 phase*/*x*x*_A*B*/*/*/*

# Generate plots based on analyzed data, but don't depend on
#   coupling and offset plots, since neat systems don't get those.
.PHONY: plot plot-energy plot-offsets plot-interface plot-hist-field plot-hist-radial plot-hist-residual plot-hist-angles
plot: plot-energy plot-offsets plot-interface plot-hist-field plot-hist-radial plot-hist-residual plot-hist-angles
plot-energy: $(PLOTS_ENERGY)
%.energy-total.png: %.energy.npy.bz2
	$(PYCULGI) plot.py $< total save
%.energy-field.png: %.energy.npy.bz2
	$(PYCULGI) plot.py $< field save
%.energy-coupl.png: %.energy.npy.bz2
	$(PYCULGI) plot.py $< coupl save
%.offsets.png: %.energy.npy.bz2
	$(PYCULGI) plot.py $< offsets save
plot-offsets: $(PLOTS_OFFSETS)
%.offsets-ang.png: %.energy.npy.bz2
	$(PYCULGI) plot.py $< offsets angles save
plot-interface: $(PLOTS_INTERFACE)
%.interface.png: %.energy.npy.bz2
	$(PYCULGI) plot.py $< interface save
plot-hist-field: $(PLOTS_HIST_FIELD)
%.hist-field-total.png: %.hist-field.npy.bz2
	$(PYCULGI) plot.py $< total save
%.hist-field-order.png: %.hist-field.npy.bz2
	$(PYCULGI) plot.py $< order save
plot-hist-radial: $(PLOTS_HIST_RADIALS)
%.hist-radial.png: %.hist-radial.npy.bz2
	$(PYCULGI) plot.py $< save
%.hist-radial.zoom.png: %.hist-radial.npy.bz2
	$(PYCULGI) plot.py $< zoom save
%.hist-radial-shell.png: %.hist-radial.npy.bz2
	$(PYCULGI) plot.py $< shell save
%.hist-radial-shell.zoom.png: %.hist-radial.npy.bz2
	$(PYCULGI) plot.py $< shell zoom save
plot-hist-residual: $(PLOTS_HIST_RES)
%.hist-residual-total.png: %.hist-residual.npy.bz2
	$(PYCULGI) plot.py $< total save
%.hist-residual-total-shell.png: %.hist-residual.npy.bz2
	$(PYCULGI) plot.py $< total shell save
%.hist-residual-order.png: %.hist-residual.npy.bz2
	$(PYCULGI) plot.py $< order save
%.hist-residual-order-shell.png: %.hist-residual.npy.bz2
	$(PYCULGI) plot.py $< order shell save
plot-hist-angles: $(PLOTS_HIST_ANGLES)
%.hist-ang.png: %.hist-ang.npy.bz2
	$(PYCULGI) plot.py $< angles save

# Generate galleries
.PHONY: gallery
gallery: gallery.html $(foreach p,$(PHASES),$(p)/gallery.html) exp/gallery.html
gallery.html: gallery.py
	$(PYCULGI) gallery.py main > gallery.html && rm -rvf Culgi.log
phase%/gallery.html: phase%/*/*/*/*.out gallery.py
	$(PYCULGI) gallery.py $* > phase$*/gallery.html && rm -rvf Culgi.log
exp/gallery.html: gallery.py exp/sem-analyzed/*/*/*.png
	$(PYCULGI) gallery.py exp > $@ && rm -rvf Culgi.log

# Generate report
.PHONY: report
report: $(NAME).pdf

# Analyze experimental data
.PHONY: exp
exp:
	$(MAKE) -C exp

# Cleanup old files (run manually)
.PHONY: cleanup
cleanup:
	"$(PYCULGI_FB)" cleanup.py

# #############################
# Remote targets (run manually)
# #############################

# Generate movies and snapshots for simulations based on Culgi outputs.
# The touching is to signify for other processes this parameter set is busy.
.PHONY: snapshots avi
snapshots: $(subst .out,.frame0001.jpg,$(SIMS_OUT_NEW))
%.frame0001.jpg: %.out %.cga
	touch $@
	-"$(PYCULGI_FB)" replay.py $< snapshots
avi: $(subst .out,.avi,$(SIMS_OUT_NEW))
%.avi: %.out %.cga
	rm -rvf $(AVIDIR)
	mkdir -p $(AVIDIR)
	touch $@
	-"$(PYCULGI_FB)" replay.py $< $(AVIDIR) save xvfb && mencoder "mf://$(AVIDIR)/*.jpg" -mf fps=10 -o $@ -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600
	rm -rvf $(AVIDIR)

# Convert Culgi output (cga/csa/ctf) to NumPy archives, and then compress
#   the targets and related files (but not the prerequisites).
# Use only the ctf file as a target, but other files also result.
.PHONY: convert
convert: $(subst .out,.ctf.npy.gz,$(SIMS_OUT_NEW))
%.ctf.npy.gz: %.out %_Inst.ctf
	touch $@
	-"$(PYCULGI_FB)" convert.py $< && gzip -f $(foreach suf,.ctf.npy .cga.npy .csa.npy,$*$(suf))

# Analyze the data (npy format, compressed). The resulting archives
#   should be much smaller than the raw data.
# Use only the energy file as a target, but other files also result.
.PHONY: analyze
analyze: $(subst .out,.energy.npy.bz2,$(SIMS_OUT))
%.energy.npy.bz2: %.out %.ctf.npy.gz
	touch $@
	-$(PYCULGI) analyze.py $< energy histograms && bzip2 -f $*.energy*.npy $*.hist*.npy
