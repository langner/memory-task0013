# ############
# Output files
# ############

# Simulation output files.
SYSTEMS = $(wildcard *x*x*_A*B*_bv?.??)
MODELS = $(foreach s,$(SYSTEMS),$(wildcard $(s)/temp*_exp*_den*_pop*))
SIMULATIONS_OUT = $(foreach m,$(MODELS),$(wildcard $(m)/*.out))
SIMULATIONS_CTF = $(foreach m,$(MODELS),$(wildcard $(m)/*.ctf))

# Analyzed data files (NumPy compressed format).
SIMULATIONS_NPZ = $(foreach m,$(MODELS),$(wildcard $(m)/*.data-analyzed.npz))

# Simulation snapshots.
SIMULATIONS_JPG = $(foreach m,$(MODELS),$(wildcard $(m)/*.jpg))

# Plots to be generated.
PLOTS_ENERGY_TOTAL = $(subst .out,.energy-total.png,$(SIMULATIONS_OUT))
PLOTS_ENERGY_FIELD = $(subst .out,.energy-field.png,$(SIMULATIONS_OUT))
PLOTS_ENERGY_COUPL = $(subst .out,.energy-coupl.png,$(SIMULATIONS_OUT))
PLOTS_HIST_RADIALS = $(subst .out,.hist-radials.png,$(SIMULATIONS_OUT))
PLOTS_ENERGY = $(PLOTS_ENERGY_TOTAL) $(PLOTS_ENERGY_FIELD) $(PLOTS_ENERGY_COUPL)
PLOTS_HISTS = $(PLOTS_HIST_RADIALS)
PLOTS_ALL = $(PLOTS_ENERGY) $(PLOTS_HISTS)

# Synchronization parameters (rsync).
SYNC_REMOTE = ~/mnt/poly/scratch/
SYNC_OPTIONS = --verbose --progress --stats --human-readable --archive --compress --update
SYNC_EXCLUDES = --exclude *.ctf --exclude *.cga --exclude *.csa --exclude *.data-raw.npz
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
.PHONY: plot plot-energy plot-hist
plot: plot-energy plot-hist
plot-energy: $(PLOTS_ENERGY)
plot-hist: $(PLOTS_HISTS)
%.energy-total.png: %.data-analyzed.npz
	python plot.py $< energy total save
%.energy-field.png: %.data-analyzed.npz
	python plot.py $< energy field save
%.energy-coupl.png: %.data-analyzed.npz
	python plot.py $< energy coupl save
%.hist-radials.png: %.data-analyzed.npz
	python plot.py $< hist radials save

# Generate the gallery if any key output files changed.
.PHONY: gallery
gallery: gallery.html
gallery.html: gallery.py $(SIMULATIONS_OUT) $(SIMULATIONS_JPG) $(PLOTS_ALL)
	python-culgi gallery.py > gallery.html

# ########################
# Remote targets (cluster)
# ########################

# Convert Culgi output (cga/csa/ctf) to NumPy compressed archives (npz).
.PHONY: convert
convert: $(subst .out,.data-raw.npz,$(SIMULATIONS_OUT))
%.data-raw.npz: %.csa %.cga %_Inst.ctf
	-"$(PYCULGI)" convert.py $(subst .csa,.out,$<)

# Analyze the data (npz format).
# The resulting archive (npz) should be much smaller.
.PHONY: analyze
analyze:  $(subst .out,.data-analyzed.npz,$(SIMULATIONS_OUT))
%.data-analyzed.npz: %.data-raw.npz
	python-culgi analyze.py $<

# Generate movies for simulations based on Culgi outputs.
.PHONY: avi
avi: $(subst .out,.avi,$(SIMULATIONS_OUT))
%.avi: %.csa %.cga
	-"$(PYCULGI)" replay.py $@ save xvfb && mencoder "mf://*.jpg" -mf fps=10 -o $@ -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600 && cp `ls *.jpg | tail -1` $(subst .avi,.jpg,$@)
	rm -rvf *.jpg


