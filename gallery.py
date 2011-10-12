# encoding: utf8

import glob
import os
import sys

from simulation import *


def controltree(sims, params):
    tree = {}
    for value in set([getattr(s,params[0]) for s in sims]):
        filtered = [s for s in sims if getattr(s,params[0]) == value]
        if params[1:]:
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
        text += "<li><b>%s:</b> %s<ul>\n" %(labels[param],value)
        if type(tree[branch]) is dict:
            text += printtree(tree[branch])
        else:
            toshow = [(printsummary(s),s) for s in tree[branch]]
            #toshow.sort()
            for summary,s in toshow:
                text += "<li>%s<ul>" %summary
                text += "\n%s\n" %printinfo(s)
                text += "</ul></li>"
        text += "</ul></li>\n"
    return text

def printsummary(sim):
    size = "%ix%ix%i" %(sim.size[0], sim.size[1], sim.size[2])
    format = "Kappa %.1f, c<sub>A</sub> %.1f, c<sub>B</sub> %.1f, exp. %.2f for %s (nchi=%.1f) in a %s box with <b>%i</b> NPs (a=%.1f)"
    params = (sim.kappa,sim.ca,sim.cb,sim.expansion,sim.polymer,sim.nchi,size,sim.population,sim.a)
    return format %params

def printinfo(sim):
    name = "%s/%s" %(sim.gallerypath,sim.outname)
    text = "<table><tr>"
    text += "<td width='300'><center>last snapshot<br/><a href='%s.jpg'><img height='250pt' src='%s.jpg' /></a></center></td>" %(name,name)
    text += "<td width='300'><center>radial distribution g(r)<a href='%s.hist-radial.png'><img border=0 height='250pt' src='%s.hist-radial.png' /></a></center></td>" %(name,name)
    text += "<td width='300'><center>radial distribution g(r)<a href='%s.hist-radial-zoom.png'><img border=0 height='250pt' src='%s.hist-radial-zoom.png' /></a></center></td>" %(name,name)
    text += "</tr></table>"
    text += "Other general: "
    text += "<a href='%s.out'>output file</a>" %name
    if sim.phase >= 4:
        text += ", <a href='%s.first.jpg'>first snapshot</a>" %name
    text += ", <a href='%s.avi'>movie (AVI)</a>" %name
    text += "<br/>"
    text += "Other energies: "
    text += "<a href='%s.energy-total.png'>total energy</a>" %name
    text += ", <a href='%s.energy-field.png'>field energy</a>" %name
    text += ", <a href='%s.energy-coupl.png'>coupling energy</a>" %name
    text += "<br/>"
    text += "Other histograms: "
    text += "<a href='%s.hist-field-total.png'>total densities</a>" %name
    text += ", <a href='%s.hist-field-order.png'>order parameters</a>" %name
    text += ", <a href='%s.hist-residual-total.png'>residual total densities</a>" %name
    text += ", <a href='%s.hist-residual-order.png'>residual order parameters</a>" %name
    text += "<br/><br/>"
    return text

control = [ "beadvolume", "temperature", "mobility", "dexcluded" ]
labels = {  "beadvolume"    : "Bead volume",
            "temperature"   : "Temperature",
            "mobility"      : "NP mobility",
            "dexcluded"     : "Demixing c<sub>A</sub>/(&rho;&kappa;+1)",
            "dcoupling"     : "Selectivity (c<sub>B</sub>-c<sub>A</sub>)"
}

phases = [ 1, 2, 3, 4, 5, 6] 
favorite = [
    # phase 1
    "phase1/64x64x1_A20B16_bv1.00/temp0.05_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0/tt11000.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0/tt11000.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0/tt11000.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0/tt11000.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0/tt11000.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0/tt11000.out",
    "phase1/64x64x1_A20B16_bv1.00/temp1.00_exp0.01_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.01_a25.0/tt11000.out",
    # phase 2
    "phase2/64x64x1_A20B16_bv1.00/temp0.01_exp1.00_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.01_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp1.00_exp0.01_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.01_a25.0/tt11000.out",
    "phase2/64x64x1_A20B16_bv1.00/temp1.00_exp0.01_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0/tt11000.out",
    # phase 3
    "phase3/64x64x1_A15B12_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.01_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.01_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.01_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.10_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.10_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.10_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca10.0_cb18.0_mob0.10_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob1.00_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob1.00_a25.0/tt11000.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi48.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    # phase 4
    "phase4/64x64x1_A15B12_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.01_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.01_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.01_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.01_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.10_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.10_a25.0/tt11000.out",
    "phase4/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.10_a25.0/tt11000.out",
    # phase 5
    "phase5/64x64x1_A15B12_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.01_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.01_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.01_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.01_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.10_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.10_a25.0/tt11000.out",
    "phase5/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.10_a25.0/tt11000.out",
    # phase 6
    "phase6/64x64x1_A15B12_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt55000.out",
    "phase6/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.01_a25.0/tt55000.out",
    "phase6/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.01_a25.0/tt55000.out",
    "phase6/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0/tt55000.out",
    "phase6/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.10_a25.0/tt55000.out",
]
favorite.sort()

HEADER = """<head>
        <script type="text/javascript" src="simpletree/simpletreemenu.js">
        </script>
        <link rel="stylesheet" type="text/css" href="simpletree/simpletree.css" />
        </head>"""

def footer(names):
    text = ""
    text += '<script type="text/javascript">\n'
    text += "//ddtreemenu.createTree(treeid, enablepersist, opt_persist_in_days (default is 1))\n"
    for name in names:
        text += 'ddtreemenu.createTree("simutree_%s", true)\n' %name
    text += "</script>\n"
    return text

description = [ "Bare NPs moving in 3D (Z&ne;0), random starting ditribution.",
                "Same as phase 1, but NPs move in 2D (Z=0)",
                "Same as phase 2, with initial density correction for high NP concentrations.",
                "Same as phase 3, with initial density equilibration.",
                "Same as phase 4, with somewhat ordered starting distribution.",
                "same as phase 5, with more ordered starting distribution, and 5x longer run time.",
]


if __name__ == "__main__":

    if "main" in sys.argv:

        print HEADER
        print "<body>"
        print "<h2>Gallery for task0013</h2>"
        print "<hr/>"
        print "<h3>Analyses of experimental images</h3>"
        print "<ul>"
        print "<li><a href='exp/gallery.html'>Separate page with analyses</a></li>"
        print "</ul"
        print "<hr/>"
        print "<h3>Simulation results</h3>"
        phases.reverse()
        for phase in phases:
            outpattern = "phase%i/*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*_ca*_cb*_mob*/t*.out" %phase
            outs = glob.glob(outpattern)
            favorite_simulations = [loadpath(out, setup=False, main=True) for out in favorite if out.count("phase%i" %phase) > 0]
            print "<h4>Phase %i (%i runs, %i favorite)</h3>" %(phase, len(outs), len(favorite_simulations))
            print "<ul>"
            print "<li>%s</li>" %description[phase-1]
            print """<li><a href="phase%i/gallery.html">Link to full gallery</a></li>""" %phase
            print "<li>Favorite runs:"
            print """<a href="javascript:ddtreemenu.flatten('simutree_phase%i', 'expand')">Expand All</a> |
                     <a href="javascript:ddtreemenu.flatten('simutree_phase%i', 'contact')">Collapse All</a>""" %(phase,phase)
            print """<ul id="simutree_phase%i" class="treeview">""" %phase
            print printtree(controltree(favorite_simulations, control))
            print "</ul></li>"
            print "</ul>"
        print footer(["phase%i" %p for p in phases])
        print "</body>"

    elif "exp" in sys.argv:

        print HEADER
        print "<body>"
        print "<h2>Gallery for Karol's analyses experimental results</h2>"
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
                    print " <a href='%s-coms.png.png'>NP centers</a>" %apath
                    print "<br/><br/>"
                    print "</ul></li>"
                print "</ul></li>"
            print "</ul></li>"
        print "</ul>"
        print footer(["sem"])
        print "</body>"

    else:

        # Phase should be passed as an integer argument.
        phase = int(sys.argv[1])

        # Get all simulations.
        outpattern = "phase%i/*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*_ca*_cb*_mob*/*.out" %phase
        outs = glob.glob(outpattern)
        outs.sort()
        simulations = [loadpath(out, setup=False) for out in outs]

        # Get favorite simulations.
        favorite_simulations = [loadpath(out, setup=False) for out in favorite if out.count("phase%i" %phase) > 0]

        print HEADER
        print "<body>"
        print "<h2>Gallery for task0013 - phase %i</h2>" %phase
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
        print footer(["favorite", "all"])
        print "</body>"
