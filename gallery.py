# encoding: utf8

import glob
import os
import sys

from simulation import *
from systems import *


def try_float(s):
    try:
        return float(s)
    except ValueError:
        return s
    except TypeError:
        return s

def controltree(sims, params):
    tree = {}
    for value in set([getattr(s,params[0]) for s in sims]):
        filtered = [s for s in sims if getattr(s,params[0]) == value]
        if params[1:] and not (params[0] == "mobility" and s.population > 0):
            tree[value,params[0]] = controltree(filtered, params[1:])
        else:
            tree[value,params[0]] = filtered
    return tree

def printtree(tree):
    text = ""
    branches = tree.keys()
    branches.sort()
    for branch in branches:
        value,param = branch
        if param == "dexcluded":
            value = "%.2f" %value
        if not (param == "mobility" and value == 0.0):
            text += "<li><b>%s:</b> %s<ul>\n" %(labels[param],value)
        else:
            text += "<li><b>Neat BCP matrix</b><ul>\n"
        if type(tree[branch]) is dict:
            text += printtree(tree[branch])
        else:
            toshow = [(printsummary(s),s) for s in tree[branch]]
            for summary,s in toshow:
                text += "<li>%s<ul>" %summary
                text += "\n%s\n" %printinfo(s)
                text += "</ul></li>"
        text += "</ul></li>\n"
    return text

def printsummary(sim):
    """
    Print the summary line that appears as the list item.
    Contents depend on whether are nanoparticles in the system.
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
    summary += " -- "
    if sim.population > 0:
        if sim.phase > 12:
            summary += "c<sub>core</sub> %.1f, " %(sim.cc)
        summary += "c<sub>A</sub> %.1f, c<sub>B</sub> %.1f" %(sim.ca, sim.cb)
    if sim.phase > 14:
        summary += ", &sigma;=%.2f" %sim.sigma

    # BCP-related parameters.
    summary += " -- bv=%.2f, nchi=%.1f, exp %.2f"%(sim.bv, sim.nchi, sim.expansion)
    if sim.phase > 9:
        summary += ", chmob=%.2f" %sim.chmob
    if sim.phase > 10:
        summary += " -- tt=%i, ts=%s" %(sim.totaltime, str(sim.timestep))

    return summary

def printinfo(sim):
    """
    Print the the contents for a single simulation run, displayed after expanding the list element.
    Perform a check on population (sim.population) for several of the things; if there are no nanoparticles,
        then the RDF cannot be displayed, for example.
    First snapshots are available only since phase 4.
    """

    name = "%s/%s" %(sim.gallerypath,sim.name)
    text = "<table><tr>"
    text += "<td width='300'><center>last snapshot<br/><a href='%s.jpg'><img height='250pt' src='%s.jpg' /></a></center></td>" %(name,name)
    if sim.population > 0:
        text += "<td width='300'><center>radial distribution g(r)<a href='%s.hist-radial.png'><img border=0 height='250pt' src='%s.hist-radial.png' /></a></center></td>" %(name,name)
        text += "<td width='300'><center>radial distribution g(r)<a href='%s.hist-radial.zoom.png'><img border=0 height='250pt' src='%s.hist-radial.zoom.png' /></a></center></td>" %(name,name)
    text += "</tr></table>"

    text += "Other visual: "
    if sim.phase >= 4:
        text += ", <a href='%s.first.jpg'>first snapshot</a>" %name
    text += ", <a href='%s.avi'>movie (AVI)</a>" %name
    text += "<br/>"

    text += "General: "
    text += "<a href='%s.out'>output file</a>" %name
    if sim.population > 0:
        text += ", <a href='%s.offsets.png'>grid offsets</a>" %name
        if sim.phase > 8:
            text += ", <a href='%s.offsets-ang.png'>angular offsets</a>" %name
            text += ", <a href='%s.hist-ang.png'>angular distribution</a>" %name
    text += "<br/>"

    text += "Energies: "
    text += "<a href='%s.energy-total.png'>total energy</a>" %name
    text += ", <a href='%s.energy-field.png'>field energy</a>" %name
    if sim.population > 0:
        text += ", <a href='%s.energy-coupl.png'>coupling energy</a>" %name
    text += "<br/>"

    text += "Field histograms: "
    text += "<a href='%s.hist-field-total.png'>total densities</a>" %name
    text += ", <a href='%s.hist-field-order.png'>order parameters</a>" %name
    text += "<br/>"

    if sim.population > 0:
        text += "Residual density histograms: "
        text += "<a href='%s.hist-residual-total.png'>residual total dens.</a>" %name
        text += ", <a href='%s.hist-residual-order.png'>residual order param.</a>" %name
        if phase > 7:
            text += ", <a href='%s.hist-residual-total-shell.png'>res. total dens. (shell)</a>" %name
            text += ", <a href='%s.hist-residual-order-shell.png'>res. order param. (shell)</a>" %name            
    text += "<br/>"

    if phase > 7:
        text += "Other g(r) plots: "
        text += "<a href='%s.hist-radial-shell.png'>g(r) for shell beads</a>" %name
        text += ", <a href='%s.hist-radial-shell.zoom.png'>g(r) for shell beads (zoomed)</a>" %name
        text += "<br/>"

    text += "<br/>"

    return text

control = [ "npname", "polymer", "kappa", "temperature", "mobility" ]
labels = {  "npname"        : "NP model",
            "polymer"       : "BCP model",
            "kappa"         : "Compressibility",
            "temperature"   : "Temperature",
            "mobility"      : "NP mobility",
}


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

def printlist(items):
    lines = "<ul>"
    for it in items:
        lines += "<li>%s</li>" %it
    lines += "</ul>"
    return lines

def footer(names, above=False):
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


if __name__ == "__main__":

    if "main" in sys.argv:

        # Header and title
        print HEADER
        print "<body>"
        print "<h2>Gallery for task0013</h2>"
        print "<hr/>"

        # Link to gallery with experimental analyses
        print "<h3>Analyses of experimental images</h3>"
        print printlist(["<a href='exp/gallery.html'>Separate page with analyses</a></li>"])
        print "<hr/>"

        # List of simulations results, with links to full galleries
        print "<h3>Simulation results</h3>"
        phases.reverse()
        print """<ul id="simutree_gallery" class="treeview">"""
        for phase in phases:
            print "<li>"
            outpattern = "phase%i/*x*x*_A*B*_bv*/temp?.??_exp?.??_den?.?_pop*/k*_nchi*/t*.out" %phase
            outs = glob.glob(outpattern)
            favorite_simulations = [loadpath(out, setup=False, main=True) for out in favorite if out.count(r"phase%i/" %phase) > 0]
            print "<h4>Phase %i (%i runs, %i favorite)</h3>" %(phase, len(outs), len(favorite_simulations))
            print printlist([
                    "%s" %description[phase-1],
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

    elif "exp" in sys.argv:

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

    else:

        # Phase should be passed as an integer argument.
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
        print "<h3>Description: %s</h3>" %description[phase-1]
        print "<h4>Favorite runs</h4>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_favorite', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_favorite', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_favorite" class="treeview">"""
        print printtree(controltree(favorite_simulations, control))
        print "</ul>"
        print "<h4>All runs</h4>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_all', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_all', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_all" class="treeview">"""
        print printtree(controltree(simulations, control))
        print "</ul>"
        print footer(["favorite", "all"], above=True)
        print "</body>"
