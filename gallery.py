# encoding: utf8

"""Script to produce a gallery of results for task0013.

These function are all quite specific to this task at hand.
"""


import glob
import os
import sys

from systems import *
from task0013 import *


def try_float(s):
    """Try to convert a float into a string."""

    try:
        return float(s)
    except ValueError:
        return s
    except TypeError:
        return s


def controltree(sims, params):
    """Build a tree from a list of simulations using params to control the structure."""

    # Start with an empty tree.
    tree = {}

    # Loop over values of the first parameter in the control list.
    for value in set([getattr(s,params[0]) for s in sims]):

        # Get all simulation where this parameter has this value.
        filtered = [s for s in sims if getattr(s,params[0]) == value]

        # If there are still parameters in the list, call this method recursively with the
        #   filtered simulations and shorter control list as input, unless the system is neat
        #   and the we've come to a NP-specific parameter.
        # Otherwise, set the tree node to the filtered list.
        if params[1:] and not (params[0] in ("mobility", "temperature") and s.population == 0):
            tree[value,params[0]] = controltree(filtered, params[1:])
        else:
            tree[value,params[0]] = filtered

    return tree


def printtree(tree, is_neat=False):
    """Build a nested unordered HTML list from a tree of simulation results.

    The flag is_neat is used to track simulations for neat systems.
    """

    text = ""

    # Get the branches and sort them.
    branches = tree.keys()
    branches.sort()

    # Loop over branches.
    for branch in branches:

        # For some reason, I chose the value to come first.
        value,param = branch

        # Check if the NP model name is up and is it's empty, and then
        #   check if the system is neat and if the parameter is NP-specific.
        is_np_none = (param == "npname") and (value == None)
        is_neat = is_neat or is_np_none
        is_neat_specific = is_neat and (param in ("mobility", "temperature"))

        # Print the branch title.
        if is_np_none:
            text += "<li><b>Neat BCP matrix</b><ul>\n"
        elif not is_neat_specific:
            text += "<li><b>%s:</b> %s<ul>\n" %(labels[param],value)            

        # If there are branches nested in this branch, recurse, otherwise
        #   append the HTML items for each simulation in this node.
        if type(tree[branch]) is dict:
            text += printtree(tree[branch], is_neat=is_neat)
        else:
            toshow = [(printsummary(s),s) for s in tree[branch]]
            for summary,s in toshow:
                text += "<li>%s<ul>" %summary
                text += "\n%s\n" %printinfo(s)
                text += "</ul></li>"

        # Add trailing closing tags if appropriate.
        if not is_neat_specific:
            text += "</ul></li>\n"

    return text


def printsummary(sim):
    """Print the summary line of a simulation.

    This line appears in the gallery tree as the list item.
    Contents depend on whether there are nanoparticles in the system.
    """

    # Population and box size.
    if sim.population > 1:
        summary = "<b>%i</b> NPs" %sim.population
    elif sim.population == 1:
        summary = "1 NP"
    else:
        summary = "Neat BCP"
    summary += " in %ix%ix%i" %(sim.size[0], sim.size[1], sim.size[2])
    if sim.population > 1:
        summary += " (a=%.1f)" %sim.a

    # NP-related parameters.
    if sim.population > 0:
        summary += " -- "
        if sim.phase > 12:
            summary += "c<sub>core</sub> %.1f, " %(sim.cc)
        summary += "c<sub>A</sub> %.1f, c<sub>B</sub> %.1f" %(sim.ca, sim.cb)
        if sim.phase > 14:
            summary += ", &sigma;=%.2f" %sim.sigma
        if sim.phase > 19:
            summary += ", stencil=%i" %sim.stencilsize

    # BCP-related parameters.
    summary += " -- bv=%.2f, nchi=%.1f, exp %.2f"%(sim.bv, sim.nchi, sim.expansion)
    if sim.phase > 9:
        summary += ", chmob=%.2f" %sim.chmob
    if sim.phase > 10:
        summary += " -- tt=%i, ts=%s" %(sim.totaltime, str(sim.timestep))

    # Population-related parameters.
    if sim.phase > 17:
        summary += " -- disp=%.2f" %(sim.disp)
    if sim.phase > 18:
        summary += ", poly=%.2f" %(sim.poly)

    return summary

def printinfo(sim):
    """
    Print the the contents for a single simulation run, displayed after expanding the list element.
    Perform a check on population (sim.population) for several of the things; if there are no nanoparticles,
        then the RDF cannot be displayed, for example.
    First snapshots are available only since phase 4, for example, and after phase 14 more snapshots are available.
    """

    # This is the upper part with images, including snapshots and the zoomed RDF.
    # Note that the name of the first and last snapshots were different before phase 14.
    # Here we assume there are always 1101 snapshots in the simulations, and that the saved
    #   snapshots are the same ones, although this might easily change in the future.
    name = "%s/%s" %(sim.gallerypath,sim.name)
    text = "<table><tr>"
    if sim.phase > 13:
        first = '%s.frame0001.jpg' %name
        last = '%s.frame1101.jpg' %name
    else:
        first = '%s.first.jpg' %name
        last = '%s.jpg' %name
    if sim.phase >= 4:
        text += "<td width='300'><center>first snapshot<br/><a href='%s'><img height='250pt' src='%s' /></a></center></td>" %(first,first)
    text += "<td width='300'><center>last snapshot<br/><a href='%s'><img height='250pt' src='%s' /></a></center></td>" %(last,last)
    if sim.population > 0:
        text += "<td width='300'><center>radial distribution g(r)<a href='%s.hist-radial.zoom.png'><img border=0 height='250pt' src='%s.hist-radial.zoom.png' /></a></center></td>" %(name,name)
    text += "</tr></table>"

    # Other snapshots and the animation.
    # After phase 16, the indexes of snapshots are different.
    text += "Other visual: "
    if sim.phase > 16:
        text += "<a href='%s.frame0011.jpg'>snapshot 11</a>, " %name
        text += "<a href='%s.frame0021.jpg'>snapshot 21</a>, " %name
        text += "<a href='%s.frame0101.jpg'>snapshot 101</a>, " %name
        text += "<a href='%s.frame0301.jpg'>snapshot 301</a>, " %name
        text += "<a href='%s.frame0501.jpg'>snapshot 501</a>, " %name
    elif sim.phase > 13:
        text += "<a href='%s.frame0045.jpg'>snapshot 45</a>, " %name
        text += "<a href='%s.frame0321.jpg'>snapshot 321</a>, " %name
    text += "<a href='%s.avi'>movie (AVI)</a>" %name
    text += "<br/>"

    # Other general stuff, including the output file and offsets.
    text += "General: "
    text += "<a href='%s.out'>output file</a>" %name
    if sim.population > 0:
        text += ", <a href='%s.offsets.png'>grid offsets</a>" %name
        if sim.phase > 8:
            text += ", <a href='%s.offsets-ang.png'>angular offsets</a>" %name
            text += ", <a href='%s.hist-ang.png'>angular distribution</a>" %name
        if sim.phase > 16:
            text += ", <a href='%s.interface.png'>NPs at interface</a>" %name
    text += "<br/>"

    # Energy plots.
    text += "Energies: "
    text += "<a href='%s.energy-total.png'>total energy</a>" %name
    text += ", <a href='%s.energy-field.png'>field energy</a>" %name
    if sim.population > 0:
        text += ", <a href='%s.energy-coupl.png'>coupling energy</a>" %name
    text += "<br/>"

    # Field histograms.
    text += "Field histograms: "
    text += "<a href='%s.hist-field-total.png'>total densities</a>" %name
    text += ", <a href='%s.hist-field-order.png'>order parameters</a>" %name
    text += "<br/>"

    # Histograms of residual densities.
    if sim.population > 0:
        text += "Residual density histograms: "
        text += "<a href='%s.hist-residual-total.png'>residual total dens.</a>" %name
        text += ", <a href='%s.hist-residual-order.png'>residual order param.</a>" %name
        if phase > 7:
            text += ", <a href='%s.hist-residual-total-shell.png'>res. total dens. (shell)</a>" %name
            text += ", <a href='%s.hist-residual-order-shell.png'>res. order param. (shell)</a>" %name            
    text += "<br/>"

    # Other RDF plots.
    if sim.pop > 0:
        text += "Other g(r) plots: "
        text += "<a href='%s.hist-radial.png'>full g(r) for cores</a>" %name
        if phase > 7:
            text += ", <a href='%s.hist-radial-shell.png'>g(r) for shell beads</a>" %name
            text += ", <a href='%s.hist-radial-shell.zoom.png'>g(r) for shell beads (zoomed)</a>" %name
            text += "<br/>"
        text += "<br/>"

    return text


def printlist(items):
    """Print a list of items."""

    return "<ul>" + "".join(["<li>%s</li>" %it for it in items]) + "</ul>"


def footer(names, above=False):
    """Print the footers necessary for the tree javascript to work."""

    text = ""
    text += '<script type="text/javascript">\n'
    text += "//ddtreemenu.createTree(treeid, enablepersist, opt_persist_in_days (default is 1))\n"
    for name in names:
        text += 'ddtreemenu.createTree("simutree_%s", true)\n' %name
    if above:
        text += 'ddtreemenu.closefolder = "../"+ddtreemenu.closefolder\n'
        text += 'ddtreemenu.openfolder = "../"+ddtreemenu.openfolder\n'
    text += "</script>\n"

    return text


# Control parameters to use.
control = [ "npname", "polymer", "kappa", "mobility", "temperature" ]

# Labels of the control parameters.
labels = {  "npname"        : "NP model",
            "polymer"       : "BCP model",
            "kappa"         : "Compressibility",
            "mobility"      : "NP mobility",
            "temperature"   : "Temperature",
}

# HTML header for the simpletree javascript code.
HEADER = """<head>
        <script type="text/javascript" src="simpletree/simpletreemenu.js">
        </script>
        <link rel="stylesheet" type="text/css" href="simpletree/simpletree.css" />
        </head>"""
HEADERabove = """<head>
        <script type="text/javascript" src="../simpletree/simpletreemenu.js">
        </script>
        <link rel="stylesheet" type="text/css" href="../simpletree/simpletree.css" />
        </head>"""


if __name__ == "__main__":

    # This refers the gallery in the root task0013 directory.
    if "main" in sys.argv:

        # We want phases in reverse order.
        phases.reverse()

        # Header and title.
        print HEADER
        print "<body><h2>Gallery for task0013</h2><hr/>"

        # Link to gallery with experimental analyses.
        print "<h3>Analyses of experimental images</h3>"
        print printlist(["<a href='exp/gallery.html'>Separate page with analyses</a></li>"])
        print "<hr/>"

        # List of simulations results, with links to full galleries.
        print "<h3>Simulation results</h3>"
        print """<ul id="simutree_gallery" class="treeview">"""

        # Loop over phases.
        for phase in phases:

            # Get all simulation output file paths for this phase.
            outpattern = "phase%i/*x*x*_A*B*_bv*/temp?.??_exp?.??_den?.?_pop*/k*_nchi*/t*.out" %phase
            outs = glob.glob(outpattern)

            # Load all simulations for this phase marked as favorite.
            favorite_simulations = [loadpath(out, setup=False, main=True) for out in favorite if out.count(r"phase%i/" %phase) > 0]

            # Print summary of phase, link to full gallery and tree of favorites.
            print "<li>"
            print "<h4>Phase %i (%i runs, %i favorite)</h3>" %(phase, len(outs), len(favorite_simulations))
            print printlist([
                    "%s" %descriptions[phase-1],
                    """<a href="phase%i/gallery.html">Full gallery on separate page</a>""" %phase,
                    "Favorite runs:<br/>"
                    + """<a href="javascript:ddtreemenu.flatten('simutree_phase%i', 'expand')">Expand All</a> |""" %phase
                    + """<a href="javascript:ddtreemenu.flatten('simutree_phase%i', 'contact')">Collapse All</a>""" %phase
                    + """<ul id="simutree_phase%i" class="treeview">""" %phase
                    + printtree(controltree(favorite_simulations, control))
                    + "</ul>"
            ])
            print "</li>"
        print "</ul>"

        # Footer
        print "<hr/>"
        print footer(["gallery"])
        print footer(["phase%i" %p for p in phases])
        print "</body>"

    # This refers to the gallery in the folder with experimental images.
    elif "exp" in sys.argv:

        # Description, link to some folders and other header-like things.
        print HEADERabove
        print "<body>"
        print "<h2>Gallery for Karol's analyses of experimental results</h2>"
        print "<h4>SEM images</h4>"
        print "<ul>"
        print "<li><a href='sem/'>Directory with all images I recieved</a>; <b>note</b>: I appended the length of the scale bar to the filename of each image</li>"
        print "<li><a href='sem-analyzed/'>Directory with all analyzed images</a></li>"
        print "<li><a href='sem-misc/'>Directory with other images (not analyzed)</a></li>"
        print "</ul>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_sem', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_sem', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_sem" class="treeview">"""

        # Benchmarks for the analysis.
        bsource = open("exp/benchmark/source").read()
        bname = "benchmark/F20-atom_res.preview"
        print "<li>Benchmark image (from <a href='%s'>%s</a>)<ul>" %(bsource,bsource)
        print "<table><tr>"
        print "<td><center>original<br/><a href='%s.png'><img src='%s.png' height='250'></a></td>" %(bname,bname)
        print "<td><center>nanoparticles<br/><a href='%s-nps.png'><img src='%s-nps.png' height='250'></a></td>" %(bname,bname)
        print "<td><center>radial distribution<br/><a href='%s-rdf.png'><img src='%s-rdf.png' height='250'></a></td>" %(bname,bname)
        print "</tr></table>"
        print "Other images:"
        print " <a href='%s-filter.png'>filtered</a>" %bname
        print " <a href='%s-thresh.png'>threshold</a>" %bname
        print " <a href='%s-seeds.png'>regional maxima</a>" %bname
        print " <a href='%s-dist.png'>distances</a>" %bname
        print " <a href='%s-coms.png'>NP centers</a>" %bname
        print "<br/><br/>"
        print "</ul></li>"

        # Now the SEM images proper that we want to present.
        pegs = [p.split('/')[-1] for p in glob.glob("exp/sem-analyzed/*k")]
        for peg in pegs:
            print "<li>%s PEG length<ul>" %peg
            Cs = [p.split('/')[-1] for p in glob.glob("exp/sem-analyzed/%s/*" %peg)]
            Cs.sort()
            Cs.append(Cs.pop(Cs.index('C10')))
            for C in Cs:
                print "<li>concentration %s<ul>" %C
                rdfs = glob.glob("exp/sem-analyzed/%s/%s/*-rdf.png" %(peg,C))
                rdfs.sort()
                for rdf in rdfs:
                    name = os.path.splitext(os.path.split(rdf)[-1])[0][:-4]
                    ipath = "sem/%s/%s/%s" %(peg,C,name)
                    apath = "sem-analyzed/%s/%s/%s" %(peg,C,name)
                    print "<li>Image %s<ul>" %name
                    print "<table><tr>"
                    print "<td><center>original<br/><a href='%s.jpg'><img src='%s.jpg' height='250'></a></td>" %(ipath,ipath)
                    print "<td><center>nanoparticles<br/><a href='%s-nps.png'><img src='%s-nps.png' height='250'></a></td>" %(apath,apath)
                    print "<td><center>radial distribution<br/><a href='%s-rdf.png'><img src='%s-rdf.png' height='250'></a></td>" %(apath,apath)
                    print "</tr></table>"
                    print "Other images:"
                    print " <a href='%s-filtered.png'>filtered</a>" %apath
                    print " <a href='%s-threshold.png'>threshold</a>" %apath
                    print " <a href='%s-coms.png'>NP centers</a>" %apath
                    print "<br/><br/>"
                    print "</ul></li>"
                print "</ul></li>"
            print "</ul></li>"
        print "</ul>"
        print footer(["sem"], above=True)
        print "</body>"

    # This finally refers to any gallery for a specific simulation phase.
    else:

        # Phase should be passed as an integer argument in this case.
        phase = int(sys.argv[1])

        # Get all simulations.
        outpattern = "phase%i/*x*x*_A*B*_bv*/temp*_exp?.??_den?.?_pop*/k*_nchi*/*.out" %phase
        outs = glob.glob(outpattern)
        outs.sort()
        simulations = [loadpath(out, setup=False) for out in outs]

        # Get favorite simulations.
        favorite_simulations = [loadpath(out, setup=False) for out in favorite if out.count("phase%i" %phase) > 0]

        print HEADERabove
        print "<body>"
        print "<h2>Gallery for task0013 - phase %i</h2>" %phase
        print "<h3>Description: %s</h3>" %descriptions[phase-1]

        # First print a tree of favorite runs.
        print "<h4>Favorite runs</h4>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_favorite', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_favorite', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_favorite" class="treeview">"""
        print printtree(controltree(favorite_simulations, control))
        print "</ul>"

        # Now print the tree with all runs.
        print "<h4>All runs</h4>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_all', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_all', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_all" class="treeview">"""
        print printtree(controltree(simulations, control))
        print "</ul>"

        print footer(["favorite", "all"], above=True)
        print "</body>"
