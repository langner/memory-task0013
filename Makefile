# ############
# Output files
# ############

# Simulation output files.
SYSTEMS = $(wildcard *x*x*_A*B*_bv?.??)
MODELS = $(foreach s,$(SYSTEMS),$(wildcard $(s)/temp*_exp*_den*_pop*))
SIMULATIONS_OUT = $(foreach m,$(MODELS),$(wildcard $(m)/*.out))
SIMULATIONS_JPG = $(foreach m,$(MODELS),$(wildcard $(m)/*.jpg))
SIMULATIONS_AVI = $(foreach m,$(MODELS),$(wildcard $(m)/*.avi))
SIMULATIONS_ENERGY = $(foreach m,$(MODELS),$(wildcard $(m)/*.energy.npy))
SIMULATIONS_HIST_FIELD = $(foreach m,$(MODELS),$(wildcard $(m)/*.hist-field.npy))
SIMULATIONS_HIST_RADIAL = $(foreach m,$(MODELS),$(wildcard $(m)/*.hist-radial.npy))
SIMULATIONS_HIST_RESIDUAL = $(foreach m,$(MODELS),$(wildcard $(m)/*.hist-residual.npy))

# Plots to be generated.
PLOTS_ENERGY_TOTAL = $(subst .energy.npy,.energy-total.png,$(SIMULATIONS_ENERGY))
PLOTS_ENERGY_FIELD = $(subst .energy.npy,.energy-field.png,$(SIMULATIONS_ENERGY))
PLOTS_ENERGY_COUPL = $(subst .energy.npy,.energy-coupl.png,$(SIMULATIONS_ENERGY))
PLOTS_HIST_FIELD_TOTAL = $(subst .hist-field.npy,.hist-field-total.png,$(SIMULATIONS_HIST_FIELD))
PLOTS_HIST_FIELD_ORDER = $(subst .hist-field.npy,.hist-field-order.png,$(SIMULATIONS_HIST_FIELD))
PLOTS_HIST_RADIAL = $(subst .hist-radial.npy,.hist-radial.png,$(SIMULATIONS_HIST_RADIAL))
PLOTS_HIST_RESIDUAL_TOTAL = $(subst .hist-residual.npy,.hist-residual-total.png,$(SIMULATIONS_HIST_RESIDUAL))
PLOTS_HIST_RESIDUAL_ORDER = $(subst .hist-residual.npy,.hist-residual-order.png,$(SIMULATIONS_HIST_RESIDUAL))
PLOTS_ENERGY = $(PLOTS_ENERGY_TOTAL) $(PLOTS_ENERGY_FIELD) $(PLOTS_ENERGY_COUPL)
PLOTS_HIST_FIELD = $(PLOTS_HIST_FIELD_TOTAL) $(PLOTS_HIST_FIELD_ORDER)
PLOTS_HIST_RESIDUAL = $(PLOTS_HIST_RESIDUAL_TOTAL) $(PLOTS_HIST_RESIDUAL_ORDER)
PLOTS_ALL = $(PLOTS_ENERGY) $(PLOTS_FIELD) $(PLOTS_RADIALS) $(PLOTS_RESIDUAL)

# Synchronization parameters (rsync).
SYNC_REMOTE = poly:scratch/
SYNC_OPTIONS = --verbose --progress --stats --human-readable --archive --compress --update
SYNC_EXCLUDES = --exclude *.ctf* --exclude *.cga* --exclude *.csa*
SYNC_DATA = /home/kml/data/

# Frame buffer parameters (xvfb).
XVFBOPTS = "-screen 0 1280x1024x24"
PYCULGI = "xvfb-run -s $(XVFBOPTS) python-culgi"

# #############
# Local targets
# #############

# Normally do plots and gallery, perhaps copy.
.PHONY: default all
default: plot gallery
all: copy plot gallery

# Data is analyzed on cluster, so synchronize local files with that.
.PHONY: copy
copy:
	rsync $(SYNC_OPTIONS) $(SYNC_EXCLUDES) $(SYNC_REMOTE)$(subst $(SYNC_DATA),,$(CURDIR))/*x*x*_A*B* .
	chmod 755 *x*x*_A*B*
	chmod 755 *x*x*_A*B*/*
	chmod 644 *x*x*_A*B*/*/*

# Generate the plots based on analyzed data.
.PHONY: plot plot-energy plot-hist-field plot-hist-radial plot-hist-residual
plot: plot-energy plot-hist-field plot-hist-radial plot-hist-residual
plot-energy: $(PLOTS_ENERGY)
%.energy-total.png %.energy-field %.energy-coupl: %.energy.npy
	python plot.py $< total save
	python plot.py $< field save
	python plot.py $< coupl save
plot-hist-field: $(PLOTS_HIST_FIELD)
%.hist-field-total.png %.hist-field-order.png: %.hist-field.npy
	python plot.py $< total save
	python plot.py $< order save
plot-hist-radial: $(PLOTS_HIST_RADIAL)
%.hist-radial.png: %.hist-radial.npy
	python plot.py $< save
plot-hist-residual: $(PLOTS_HIST_RESIDUAL)
%.hist-residual-total.png %.hist-residual-order.png: %.hist-residual.npy
	python plot.py $< total save
	python plot.py $< order save

# Generate the gallery if any key output files changed.
.PHONY: gallery
gallery: gallery.html
gallery.html: gallery.py $(SIMULATIONS_OUT) $(SIMULATIONS_JPG) $(SIMULATIONS_AVI) $(PLOTS_ALL)
	python-culgi gallery.py > gallery.html

# ########################
# Remote targets (cluster)
# ########################

# Convert Culgi output (cga/csa/ctf) to NumPy compressed archives (npz).
.PHONY: convert
convert: $(subst .out,.ctf.npy,$(SIMULATIONS_OUT)) $(subst .out,.csa.npy,$(SIMULATIONS_OUT)) $(subst .out,.cga.npy,$(SIMULATIONS_OUT))
%.ctf.npy %.csa.npy %.cga.npy: %.out %_Inst.ctf %.csa %.cga
	-"$(PYCULGI)" convert.py $<

# Analyze the data (npz format).
# The resulting archive (npz) should be much smaller than the raw data.
.PHONY: analyze
analyze:  $(subst .out,.energy.npy,$(SIMULATIONS_OUT)) $(subst .out,.hist-field.npy,$(SIMULATIONS_OUT)) $(subst .out,.hist-radial.npy,$(SIMULATIONS_OUT)) $(subst .out,.hist-residual.npy,$(SIMULATIONS_OUT))
%.energy.npy %.hist-field.npy %.hist-radial.npy %.hist-residual.npy: %.out %.ctf.npy %.csa.npy %.cga.npy
	-python-culgi analyze.py $<

# Generate movies and snapshots for simulations based on Culgi outputs.
.PHONY: avi
avi: $(subst .out,.avi,$(SIMULATIONS_OUT))
%.avi: %.csa %.cga
	-"$(PYCULGI)" replay.py $@ save xvfb && mencoder "mf://*.jpg" -mf fps=10 -o $@ -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600 && cp `ls *.jpg | tail -1` $(subst .avi,.jpg,$@)
	rm -rvf *.jpg


