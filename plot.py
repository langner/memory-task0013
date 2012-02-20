import bz2
import sys

import numpy
import pylab

from analyze import PhaseFrames


def correction(X,Y,R):
    ang1 = 2 * numpy.pi
    ang2 = 2 * numpy.pi - 2.0
    ang3 = 2 * numpy.pi - 2.0
    ang4 = 3 * numpy.pi / 2.0 - 2.0
    V1 = (0.5*X - R)*(0.5*Y - R)
    V2 = R*(0.5*X - R)
    V3 = R*(0.5*Y - R)
    V4 = R**2
    ref = (V1+V2+V3+V4) * 2 * numpy.pi
    return (ang1*V1 + ang2*V2 + ang3*V3 + ang4*V4) / ref


if __name__ == "__main__":

    fpath = sys.argv[1].strip()
    phase,model,system,params,run = fpath.split('/')
    phase = int(phase[5:])
    size = map(int,model.split("_")[0].split('x'))
    population = int(system.split('_')[-1][3:])
    if phase == 8 and "shell" in sys.argv:
        population *= 8
    data = numpy.load(bz2.BZ2File(fpath))

    if "energy.npy.bz2" in fpath:

        froot = fpath.replace(".energy.npy.bz2","")
        E = data

        if "total" in sys.argv:
            instant = data[:,1] + data[:,3] + data[:,4]  + data[:,5]  + data[:,6]
            average = data[:,7] + data[:,9] + data[:,10] + data[:,11] + data[:,12]
            datas = [ "instant.", "average" ]
            dataset = { "instant.": instant, "average": average }
            plotfname = "%s.energy-total.png" %froot
            ylabel = "total free energy"
        elif "field" in sys.argv:
            instant = data[:,3] + data[:,4]  + data[:,5]
            average = data[:,9] + data[:,10] + data[:,11]
            datas = [ "instant.", "average" ]
            dataset = { "instant.": instant, "average": average }
            plotfname = "%s.energy-field.png" %froot
            ylabel = "field free energy"
        elif "coupl" in sys.argv and population > 0:
            instant = data[:,6]
            average = data[:,12]
            datas = [ "instant.", "average" ]
            dataset = { "instant.": instant, "average": average }
            plotfname = "%s.energy-coupl.png" %froot
            ylabel = "coupling free energy"
        elif "offsets" in sys.argv and phase >= 7 and population > 0:
            if phase == 7:
                mean = data[:,19]
                var = data[:,20]
                skew = data[:,21]
                kurtosis = data[:,22]
                datas = [ "mean", "variance", "skewness", "kurtosis" ]
                dataset = { "mean": mean, "variance": var, "skewness": skew, "kurtosis": kurtosis }
            else:
                mean = data[:,19]
                mean_shell = data[:,20]
                var = data[:,21]
                var_shell = data[:,22]
                skew = data[:,23]
                skew_shell = data[:,24]
                kurtosis = data[:,25]
                kurtosis_shell = data[:,26]
                datas = [ "mean", "mean (shell)", "variance", "variance (shell)", "skewness", "skewness (shell)", "kurtosis", "kurtosis (shell)" ]
                dataset = { "mean": mean, "mean (shell)": mean_shell, "variance": var, "variance (shell)": var_shell,
                            "skewness": skew, "skewness (shell)": skew_shell, "kurtosis": kurtosis, "kurtosis (shell)": kurtosis_shell }
            plotfname = "%s.offsets.png"  %froot
            ylabel = "offset relative to grid midpoint"
        else:
            sys.exit(0)

        for ds in datas:
            pylab.plot(data[:,0], dataset[ds], label=ds)

        ymin = min([min(dataset[ds][100:]) for ds in dataset])
        ymax = max([max(dataset[ds][100:]) for ds in dataset])
        dy = ymax - ymin
        if "field" in sys.argv:
            ymin = ymax - dy/1.05
            ymax = ymax + dy/100
        else:
            ymin = ymin - dy/100
            ymax = ymin + dy/1.05
        if "offsets" in sys.argv:
            ymin = -1.3
            ymax = 0.6

        pylab.xlabel("time step")
        pylab.ylabel(ylabel)

        xmin = data[0,0]
        xmax = data[-1,0]*1.1
        pylab.xlim(xmin,xmax)
        pylab.ylim(ymin,ymax)
        pylab.legend(loc="best")

        if "offsets" in sys.argv:
            pylab.axhline(y=0.5, xmin=data[0,0], xmax=data[-1,0], linestyle='--', color='gray')
            pylab.axhline(y=1.0/12, xmin=data[0,0], xmax=data[-1,0], linestyle='--', color='gray')
            pylab.axhline(y=0.0, xmin=data[0,0], xmax=data[-1,0], linestyle='--', color='gray')
            pylab.axhline(y=-1.2, xmin=data[0,0], xmax=data[-1,0], linestyle='--', color='gray')
            pylab.text(xmax*1.01, 0.5, "0.5", fontsize=12, horizontalalignment="left", verticalalignment="center")
            pylab.text(xmax*1.01, 1.0/12, "1/12", fontsize=12, horizontalalignment="left", verticalalignment="center")
            pylab.text(xmax*1.01, 0.0, "0.0", fontsize=12, horizontalalignment="left", verticalalignment="center")
            pylab.text(xmax*1.01, -1.2, "-1.2", fontsize=12, horizontalalignment="left", verticalalignment="center")

    if "hist-field.npy.bz2" in fpath:
        froot = fpath.replace(".hist-field.npy.bz2","")

        if "total" in sys.argv:
            plotfname = froot+".hist-field-total.png"
            xlabel = "overall field density"
            xmin, xmax = 0.0, 1.5
            data = data[:,0]

        if "order" in sys.argv:
            plotfname = froot+".hist-field-order.png"
            xlabel = "overall order parameter"
            xmin, xmax = -1.5, 1.5
            data = data[:,1]

    if "hist-radial.npy.bz2" in fpath:

        froot = fpath.replace(".hist-radial.npy.bz2","")

        if phase > 7:
            if "shell" in sys.argv:
                data = data[:,1,:]
                plotfname = froot+".hist-radial-shell"
            else:
                data = data[:,0,:]
                plotfname = froot+".hist-radial"
        else:
            plotfname = froot+".hist-radial"

        # There were no shells before phase 8
        if phase <= 7 and "shell" in sys.argv:
            print "No shells before phase 8"
            exit(0);

        if "zoom" in sys.argv:
            plotfname += ".zoom.png"
        else:
            plotfname += ".png"

        # The X axis is fixed from analysis, but we change the plotted scale later
        xlabel = "NP-NP distance"
        xmin = 0.0
        xmax = 16.0

    if "hist-residual.npy.bz2" in fpath:

        froot = fpath.replace(".hist-residual.npy.bz2","")

        if "total" in sys.argv:
            plotfname = froot+".hist-residual-total.png"
            xlabel = "residual field density"
            xmin, xmax = 0.0, 1.5
            data = data[:,0]

        if "order" in sys.argv:
            plotfname = froot+".hist-residual-order.png"
            xlabel = "residual order parameter"
            xmin, xmax = -1.5, 1.5
            data = data[:,1]


    if "hist" in fpath:

        nsamples = PhaseFrames[phase][1]
        frames = numpy.load(bz2.BZ2File(froot+".hist-frames.npy.bz2")).tolist()
        nbins = data.shape[1]

        dx = (xmax-xmin)/nbins
        xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)

        Npairs = population*(population-1)/2
        if "hist-radial.npy" in fpath:

            # Cut out the strong peaks due to internal structure
            # Must also adjust the effective normalization factor
            # This very specific hack works for phase 8
            if (phase == 8) and ("shell" in sys.argv):
                a, b = 0.2828, 0.4
                r1 = numpy.sqrt( (b-a)**2 + a**2 )
                r2 = numpy.sqrt(2.0)*b
                r3 = numpy.sqrt( (b+a)**2 + a**2 )
                r4 = 2*b
                b1 = int(r1//dx)
                b2 = int(r2//dx)
                b3 = int(r3//dx)
                b4 = int(r4//dx)
                data[:,b1] -= population
                data[:,b2] -= population
                data[:,b3] -= population
                data[:,b4] -= population/2
                Npairs -= population*25

            strips = correction(size[0],size[1],xrange) * 2 * numpy.pi * xrange * dx
            factor = size[0] * size[1] / strips / Npairs
            hist = factor * data[frames.index(0)]
            pylab.plot(xrange, hist, label="frame 0")
        else:
            factor = 1.0 / sum(data[0])

        hist = factor * numpy.sum([data[frames.index(1+i)] for i in range(nsamples)], axis=0) / nsamples
        pylab.plot(xrange, hist, label="frames 1-%i" %(1+nsamples))
        hist = factor * numpy.sum([data[frames.index(101+i)] for i in range(nsamples)], axis=0) / nsamples
        pylab.plot(xrange, hist, label="frames 101-%i" %(101+nsamples))
        hist = factor * numpy.sum([data[frames.index(1001+i)] for i in range(nsamples)], axis=0) / nsamples
        pylab.plot(xrange, hist, label="frames 1001-%i" %(1001+nsamples))
        hist = factor * numpy.sum([data[frames.index(10001+i)] for i in range(nsamples)], axis=0) / nsamples
        pylab.plot(xrange, hist, label="frames 10001-%i" %(10001+nsamples))
        if 50001 in frames:
            hist = factor * numpy.sum([data[frames.index(50001+i)] for i in range(nsamples)], axis=0) / nsamples
            pylab.plot(xrange, hist, label="frames 50001-%i" %(50001+nsamples))
        pylab.xlabel(xlabel)
        pylab.ylabel("probability density")

        pylab.grid()
        pylab.legend()

        # We are most interested in small radial distances
        if "hist-radial.npy" in fpath:
            xmin = 0.0
            if "zoom" in sys.argv:
                if phase > 7:
                    xmin = 1.35
                    xmax = 3.7
                else:
                    xmax = 2.2
            else:
                xmax = 16.0
            pylab.xlim(xmin,xmax)

    if "save" in sys.argv:
        pylab.savefig(plotfname)
    else:
        pylab.show()