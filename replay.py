"""Replay and save visualization for task0013."""


import os
import sys

from pyculgi import *

from task0013 import gen_shell_beads, loadpath


if __name__ == "__main__":

    # First argument must be a path to load into a simulation.
    # Load the simulation, but don't do the setup yet.
    sim = loadpath(sys.argv[1], setup=False)

    # Second argument must be a temporary directory that we can use to store lots of files.
    tmpdir = sys.argv[2]

    # Some flags.
    tosave = "save" in sys.argv
    issnapshots = "snapshots" in sys.argv
    isoffscreen = "xvfb" in sys.argv or "offscreen" in sys.argv

    # Make sure that simulation has finished (if we want to save).
    # Note that from phase 4, there are two runs for each simulation.
    isfinished = "Time used" in open(sim.outpath).read().strip().split('\n')[-1]
    if (issnapshots or tosave) and not isfinished:
        print "This simulation has not finished."
        sys.exit(1)

    # Since we are serious about this run, do the setup now.
    sim.setup()

    # Tweak colors and display configuration.
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

    # Create and load graphics object.
    graphics = GraphicsManager.CreateGraphics()
    graphics.AddViewable(sim.box)
    graphics.AddViewable(sim.bcp_A)
    graphics.AddViewable(sim.bcp_B)

    # Set bead dispay radii and add them to the graphics object.
    # Remember that after phase 7 the nanoparticles are colloids.
    # After phase 18 there is a whole range of polydisperse shell beads,
    #   and we want to set the radius to zero for all of these.
    for i in GPERange(0, sim.population, 1, "<"):
        if sim.phase <= 7:
            npi = sim.box.GetSoftCoreMoleculeCmds("np", i)
            npi.SetBeadDisplayRadius("P", 0.75)
        else:
            npi = sim.box.GetSoftCoreColloidCmds("np", i)
            npi.SetBeadDisplayRadius("C", 1.5)
            npi.SetBeadDisplayRadius("S", 0.0)
            if sim.phase > 18:
                for digits,sign,delta,name in gen_shell_beads():
                    npi.SetBeadDisplayRadius(name, 0.0)
            
        graphics.AddViewable(npi)

    # Paths to cga and csa archive files.
    cga = "%s/%s.cga" %(sim.path,sim.name)
    if sim.population > 0:
        csa = "%s/%s.csa" %(sim.path,sim.name)

    # Create and configure that archive reader.
    archive = UtilitiesManager.CreateArchiveReader()
    if sim.population > 0:
        archive.SetSystem(sim.box, cga, csa)
    else:
        archive.SetSystem(sim.box, cga)
    archive.AddGraphics(graphics)

    # Final configuration actions for the graphics object.
    if isoffscreen:
        graphics.SetOffScreenRenderingOn()
    graphics.SetDisplay3DAxesOff()
    graphics.SetSize(500,500)
    graphics.SetPosition(100,100)
    graphics.Zoom(1.3)

    # Nuber of frames in the archive.
    N = archive.GetNumberOfFrames()

    # If we want snaphots, save just four analogously to experimental images.
    # Notice how the frames are scaled by 2, 14 and 48, which are the experimemntal
    #   annealing times for which the SEM images were created.
    # The little dance with temporary files is needed, because the first snapshot
    #   is touched in a makefile for logistic purposes.
    if issnapshots:

        # Which frames to save as snapshots.
        # At first this was done linearly, with frames in progression scaled linearly
        #   with the experimental times, but now we use different frames, because
        #   the correspondence is not linear, and different for NPs and the BCP.
        snaps = [1, 11, 21, 101, 301, 501, 1101]
        snaps.sort()
        for i in snaps:

            print "Saving snapshot %i..." %i
            archive.LoadFrame(i)
            fname = "%s/%s.frame%s" %(sim.path, sim.name, str(i).zfill(4))
            graphics.WriteJPEG(fname+"_tmp")
            if os.path.exists(fname+".jpg"):
                os.remove(fname+".jpg")
            os.rename(fname+"_tmp.jpg", fname+".jpg")

    # If we want to save the movie, we need to figure out the FPS and such.
    # The movie will have to parts, a slow part (10s) that uses each frame,
    #   and a second part (20s) that uses all remaining frames, sampling them
    #   with the corresponding frequency at the same FPS.
    elif tosave:

        fps = 10.0
        tslow = 10.0
        tnormal = 20.0
        ifreq = int((N-tslow*fps)/(fps*tnormal))
        iframe = 1
        inum = 1
        for i in range(int(tslow*fps)):
            print "Loading frame %i..." %iframe
            archive.LoadFrame(iframe)
            num = ("%i" %inum).zfill(4)
            graphics.SetInputText("")
            graphics.WriteJPEG("%s/frame_%s" %(tmpdir,num))
            iframe += 1
            inum += 1
        while iframe < N:
            print "Loading frame %i..." %iframe
            archive.LoadFrame(iframe)
            num = ("%i" %inum).zfill(4)
            graphics.SetInputText("")
            graphics.WriteJPEG("%s/frame_%s" %(tmpdir,num))
            iframe += ifreq
            inum += 1

    # Otherwise, just play the movie.
    else:
    
        archive.Play()
        Culgi.StartEventLoop()
