# encoding: utf8

import glob

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
    name = "%s/%s" %(sim.outpath,sim.outname)
    text = "<table><tr>"
    text += "<td width='300'><center>last snapshot<br/><a href='%s.jpg'><img height='250pt' src='%s.jpg' /></a></center></td>" %(name,name)
    text += "<td width='300'><center>radial distribution g(r)<a href='%s.hist-radials.png'><img border=0 height='250pt' src='%s.hist-radials.png' /></a></center></td>" %(name,name)
    text += "<td width='300'><center>residual order params.<a href='%s.res-orders.png'><img border=0 height='250pt' src='%s.res-orders.png' /></a></center></td>" %(name,name)
    text += "</tr></table>"
    text += "Other general: "
    text += "<a href='%s.out'>output file</a>" %name
    text += ", <a href='%s.avi'>movie (AVI)</a>" %name
    text += "<br/>"
    text += "Other energies: "
    text += "<a href='%s/%s.energy-total.png'>total energy</a>" %(sim.outpath,sim.outname)
    text += ", <a href='%s/%s.energy-field.png'>field energy</a>" %(sim.outpath,sim.outname)
    text += ", <a href='%s/%s.energy-coupl.png'>coupling energy</a>" %(sim.outpath,sim.outname)
    text += "<br/>"
    text += "Other histograms: "
    text += "<a href='%s/%s.hist-totals.png'>total densities</a>" %(sim.outpath,sim.outname)
    text += ", <a href='%s/%s.hist-orders.png'>order parameters</a>" %(sim.outpath,sim.outname)
    text += ", <a href='%s/%s.res-totals.png'>residual total densities</a>" %(sim.outpath,sim.outname)
    text += "<br/><br/>"
    return text


outpattern = "*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*_ca*_cb*_mob*.out"
outs = glob.glob(outpattern)
simulations = [loadpath(out) for out in outs]
Nsimulations = len(simulations)

control = [ "beadvolume", "temperature", "mobility", "dexcluded", "dcoupling" ]
labels = {  "beadvolume"    : "Bead volume",
            "temperature"   : "Temperature",
            "mobility"      : "NP mobility",
            "dexcluded"     : "Demixing (c<sub>A</sub>-&kappa;)",
            "dcoupling"     : "Selectivity (c<sub>B</sub>-c<sub>A</sub>)"
}
Ncontrol = len(control)

selected = [
    # temperature 0.05
    "64x64x1_A20B16_bv1.00/temp0.05_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.05_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    # temperature 0.10
    # mobility 0.10
    # selectivity 2.0
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb18.0_mob0.10_a25.0.out",
    # selectivity 8.0
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb24.0_mob0.10_a25.0.out",
    # selectivity 16.0
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb32.0_mob0.10_a25.0.out",
    # mobility 1.00
    # selectivity 2.0
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb18.0_mob1.00_a25.0.out",
    # selectivity 8.0
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop200/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop500/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop1000/k15.0_nchi24.0_ca16.0_cb24.0_mob1.00_a25.0.out",
    # selectivity 16.0
    "64x64x1_A20B16_bv1.00/temp0.10_exp0.10_den1.0_pop100/k15.0_nchi24.0_ca16.0_cb32.0_mob1.00_a25.0.out",
]
selected_simulations = [loadpath(out) for out in selected]

print """<head>
<script type="text/javascript" src="simpletree/simpletreemenu.js">
/***********************************************
* Simple Tree Menu- © Dynamic Drive DHTML code library (www.dynamicdrive.com)
* This notice MUST stay intact for legal use
* Visit Dynamic Drive at http://www.dynamicdrive.com/ for full source code
***********************************************/
</script>
<link rel="stylesheet" type="text/css" href="simpletree/simpletree.css" />
</head>"""

print "<body>"

print """<h3>Gallery of task0013a simulations</h3>"""

print """
<h4>Selected simulations</h4>
<a href="javascript:ddtreemenu.flatten('simutree_selected', 'expand')">Expand All</a> | <a href="javascript:ddtreemenu.flatten('simutree_selected', 'contact')">Collapse All</a>
<ul id="simutree_selected" class="treeview">"""

print printtree(controltree(selected_simulations, control))

print """</ul>
<br/><hr/><br/>
<h4>All simulations</h4>
<a href="javascript:ddtreemenu.flatten('simutree', 'expand')">Expand All</a> | <a href="javascript:ddtreemenu.flatten('simutree', 'contact')">Collapse All</a>
<ul id="simutree" class="treeview">"""

print printtree(controltree(simulations, control))

print "</ul>"

print """<script type="text/javascript">
//ddtreemenu.createTree(treeid, enablepersist, opt_persist_in_days (default is 1))
ddtreemenu.createTree("simutree_selected", true)
ddtreemenu.createTree("simutree", true)
</script>"""

print "</body>"
