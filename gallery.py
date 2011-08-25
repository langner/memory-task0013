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
    text = "View <a href='%s/%s.out'>simulation output</a>" %(sim.outpath,sim.outname)
    text += ", <a href='%s/%s.avi'>animation (AVI)</a><br/>" %(sim.outpath,sim.outname)
    text += "<img height='250pt' src='%s/%s.jpg' />" %(sim.outpath,sim.outname)
    text += "<img height='250pt' src='%s/%s.energy.png' />" %(sim.outpath,sim.outname)
    text += "<br/><br/>"
    return text


outpattern = "*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*_ca*_cb*_mob*.out"
outs = glob.glob(outpattern)
simulations = [loadpath(out) for out in outs]
Nsimulations = len(simulations)

control = [ "beadvolume", "temperature", "mobility", "dcoupling" ]
labels = {  "beadvolume"	: "Bead volume",
            "temperature"	: "Temperature",
            "mobility"        	: "NP mobility",
            "dcoupling"		: "Selectivity",
            "kappa"		: "Compressibility",
}
Ncontrol = len(control)

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

print """<h4>Gallery of task0013 simulations</h4>
<a href="javascript:ddtreemenu.flatten('simutree', 'expand')">Expand All</a> | <a href="javascript:ddtreemenu.flatten('simutree', 'contact')">Collapse All</a>
<ul id="simutree" class="treeview">"""

print printtree(controltree(simulations, control))

print "</ul>"

print """<script type="text/javascript">
//ddtreemenu.createTree(treeid, enablepersist, opt_persist_in_days (default is 1))
ddtreemenu.createTree("simutree", true)
</script>"""

print "</body>"
