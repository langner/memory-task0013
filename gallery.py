# encoding: utf8

import glob
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
    format = "Kappa %.1f, c<sub>A</sub> %.1f, exp. %.2f for %s (nchi=%.1f) in a %s box with <b>%i</b> NPs (a=%.1f)"
    params = (sim.kappa,sim.ca,sim.expansion,sim.polymer,sim.nchi,size,sim.population,sim.a)
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
    text += ", <a href='%s.avi'>movie (AVI)</a>" %name
    text += "<br/>"
    text += "Other energies: "
    text += "<a href='%s.energy-total.png'>total energy</a>" %name
    text += ", <a href='%senergy-field.png'>field energy</a>" %name
    text += ", <a href='%senergy-coupl.png'>coupling energy</a>" %name
    text += "<br/>"
    text += "Other histograms: "
    text += "<a href='%s.hist-field-total.png'>total densities</a>" %name
    text += ", <a href='%s.hist-field-order.png'>order parameters</a>" %name
    text += ", <a href='%s.hist-residual-total.png'>residual total densities</a>" %name
    text += ", <a href='%s.hist-residual-order.png'>residual order parameters</a>" %name
    text += "<br/><br/>"
    return text

control = [ "beadvolume", "temperature", "mobility", "dexcluded", "dcoupling" ]
labels = {  "beadvolume"    : "Bead volume",
            "temperature"   : "Temperature",
            "mobility"      : "NP mobility",
            "dexcluded"     : "Demixing (c<sub>A</sub>-&kappa;)",
            "dcoupling"     : "Selectivity (c<sub>B</sub>-c<sub>A</sub>)"
}

selected = [
    "phase1/64x64x1_A20B16_bv1.00/temp0.05_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "phase1/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb32.0_mob1.00_a25.0.out",

    # temperature 0.05
    "phase2/64x64x1_A20B16_bv1.00/temp0.05_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.05_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    # temperature 0.10

    # mobility 0.10
    # selectivity 2.0
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    # selectivity 8.0
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    # selectivity 16.0
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    # mobility 1.00
    # selectivity 2.0
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    # selectivity 8.0
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    # selectivity 16.0
    "phase2/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb32.0_mob1.00_a25.0.out",

    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.01_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.01_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.01_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.01_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob0.10_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca7.0_cb14.0_mob0.10_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob0.10_a50.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca9.0_cb16.0_mob0.10_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca10.0_cb18.0_mob0.10_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca6.0_cb12.0_mob1.00_a25.0.out",
    "phase3/64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca8.0_cb16.0_mob1.00_a25.0.out",
]

HEADER = """<head>
        <script type="text/javascript" src="simpletree/simpletreemenu.js">
        /***********************************************
        * Simple Tree Menu- © Dynamic Drive DHTML code library (www.dynamicdrive.com)
        * This notice MUST stay intact for legal use
        * Visit Dynamic Drive at http://www.dynamicdrive.com/ for full source code
        ***********************************************/
        </script>
        <link rel="stylesheet" type="text/css" href="simpletree/simpletree.css" />
        </head>"""

FOOTER = """<script type="text/javascript">
        //ddtreemenu.createTree(treeid, enablepersist, opt_persist_in_days (default is 1))
        ddtreemenu.createTree("simutree_phase1", true)
        ddtreemenu.createTree("simutree_phase2", true)
        ddtreemenu.createTree("simutree_phase3", true)
        ddtreemenu.createTree("simutree_selected", true)
        ddtreemenu.createTree("simutree_all", true)
        </script>"""

description = [ "Bare NPs, no Z constraint, random initial distribution",
                "Bare NPs, Z constraint, random initial distribution",
                "Bare NPs, Z constraint, random initial distribution"
]


if __name__ == "__main__":

    if "main" in sys.argv:

        print HEADER
        print "<body>"
        print "<h2>Gallery for task0013</h2>"
        for phase in 1,2,3:
            outpattern = "phase%i/*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*_ca*_cb*_mob*.out" %phase
            outs = glob.glob(outpattern)
            selected_simulations = [loadpath(out, setup=False, main=True) for out in selected if out.count("phase%i" %phase) > 0]
            print "<h3>Phase %i (%i runs, %i selected)</h3>" %(phase, len(outs), len(selected_simulations))
            print "<ul>"
            print "<li>Description: %s</li></br>" %description[phase-1]
            print """<li><a href="phase%i/gallery.html">Link to full gallery</a></li><br/>""" %phase
            print "<li>Selected runs:"
            print """<a href="javascript:ddtreemenu.flatten('simutree_phase%i', 'expand')">Expand All</a> |
                     <a href="javascript:ddtreemenu.flatten('simutree_phase%i', 'contact')">Collapse All</a>""" %(phase,phase)
            print """<ul id="simutree_phase%i" class="treeview">""" %phase
            print printtree(controltree(selected_simulations, control))
            print "</ul></li><br/>"
            print "</ul>"
        print FOOTER
        print "</body>"

    else:

        # Phase should be passed as an integer argument.
        phase = int(sys.argv[1])

        # Get all simulations.
        outpattern = "phase%i/*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*_ca*_cb*_mob*.out" %phase
        outs = glob.glob(outpattern)
        simulations = [loadpath(out, setup=False) for out in outs]

        # Get selected simulations.
        selected_simulations = [loadpath(out, setup=False) for out in selected if out.count("phase%i" %phase) > 0]

        print HEADER
        print "<body>"
        print "<h2>Gallery for task0013 - phase %i</h2>" %phase
        print "<h4>Selected runs</h4>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_selected', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_selected', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_selected" class="treeview">"""
        print printtree(controltree(selected_simulations, control))
        print "</ul>"
        print "<h4>All runs</h4>"
        print """<a href="javascript:ddtreemenu.flatten('simutree_all', 'expand')">Expand All</a> |
                 <a href="javascript:ddtreemenu.flatten('simutree_all', 'contact')">Collapse All</a>"""
        print """<ul id="simutree_all" class="treeview">"""
        print printtree(controltree(simulations, control))
        print "</ul>"
        print FOOTER
        print "</body>"
