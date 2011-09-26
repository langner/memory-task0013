import sys

from simulation import *


if __name__ == "__main__":

    sim = loadpath(sys.argv[1])
    out = "%s/%s.out" %(sim.outpath,sim.outname)

    if "save" in sys.argv and not "Time used" in open(out).read():

        print "This simulation has not finished."
        sys.exit(1)

    sim.bcp_A.ClearColorMap()
    sim.bcp_B.ClearColorMap()
    sim.bcp_A.AddVolumeDisplayRGBPoint(0.0, 0, 0, 256)
    sim.bcp_A.AddVolumeDisplayRGBPoint(0.5, 0, 256, 0)
    sim.bcp_A.AddVolumeDisplayRGBPoint(1.0, 256, 0, 0)
    sim.bcp_B.AddVolumeDisplayRGBPoint(0.0, 256, 0, 0)
    sim.bcp_B.AddVolumeDisplayRGBPoint(0.5, 0, 256, 0)
    sim.bcp_B.AddVolumeDisplayRGBPoint(1.0, 0, 0, 256)
    sim.bcp_A.SetRayCastingSampleDistance(0.1)
    sim.bcp_B.SetRayCastingSampleDistance(0.1)

    graphics = GraphicsManager.CreateGraphics()
    graphics.AddViewable(sim.box)
    graphics.AddViewable(sim.bcp_A)
    graphics.AddViewable(sim.bcp_B)
    for i in GPERange(0, sim.population, 1, "<"):
        npi = sim.box.GetSoftCoreMoleculeCmds("np", i)
        npi.SetBeadDisplayRadius("P", 0.75)
        graphics.AddViewable(npi)

    cga = "%s/%s.cga" %(sim.outpath,sim.outname)
    csa = "%s/%s.csa" %(sim.outpath,sim.outname)

    archive = UtilitiesManager.CreateArchiveReader()
    if "xvfb" in sys.argv or "offscreen" in sys.argv:
        graphics.SetOffScreenRenderingOn()

    archive.SetSystem(sim.box, cga, csa)
    archive.AddGraphics(graphics)

    graphics.SetDisplay3DAxesOff()
    graphics.SetInputText("")
    graphics.SetSize(500,500)
    graphics.SetPosition(100,100)
    graphics.Zoom(1.3)

    if "save" in sys.argv:

        N = archive.GetNumberOfFrames()
        fps = 10.0
        tslow = 10.0
        tnormal = 20.0
        ifreq = int((N-tslow*fps)/(fps*tnormal))
        iframe = 1
        inum = 1
        for i in range(int(tslow*fps)):
            archive.LoadFrame(iframe)
            num = ("%i" %inum).zfill(4)
            graphics.WriteJPEG("frame_%s" %num)
            iframe += 1
            inum += 1
        while iframe < N:
            archive.LoadFrame(iframe)
            num = ("%i" %inum).zfill(4)
            graphics.WriteJPEG("frame_%s" %num)
            iframe += ifreq
            inum += 1

    elif "last" in sys.argv:

        archive.LoadFrame(1)
        Culgi.StartEventLoop()

    else:
    
        archive.Play()
        Culgi.StartEventLoop()
