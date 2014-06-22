"""
Todo:
    - only average trajectory over certain time interval, which maybe can be given ...
"""

from numpy import *
from matplotlib.pyplot import * #xlabel,ylabel,plot,legend,show,step,connect,title,setp
from matplotlib import animation
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument("rundir", help="Directory which contains 'simOutput'")
parser.add_argument("-sa", "--saveani", help="save animation of 2D Simulation", action="store_true")
parser.add_argument("-e", "--initenergies", help="", action="store_true")
parser.add_argument("-s", "--stats", help="", action="store_true")
parser.add_argument("-a", "--particlemovement", help="", action="store_true")
parser.add_argument("-l", "--lastconfig", help="", action="store_true")
parser.add_argument("-d", "--density", help="", action="store_true")
parser.add_argument("-np", "--numofparticles", help="", type=int, required=True )
parser.add_argument("--tmax", help="", type=float)
args = parser.parse_args()
N_Particles = args.numofparticles

oldSimdataAvailable = os.path.isfile( str(args.rundir)+"/simdata.dat" )

# loadtext/genfromtxt always loads the file row-wise in an array data[0] being the first row

# Energy Distributions for uniformly random start positions
if args.initenergies:
    figure()
    data = genfromtxt( str(args.rundir)+"/CellEnergies.dat", dtype='float', comments='#' )
    step( data[0],data[1], where='mid', label="Coulomb" )
    step( data[0],data[2], where='mid', label="Spheres" )
    step( data[0],data[3], where='mid', label="CIC" )
    legend()

# E,T,L,V,P for random initial positions going into equilibrium through stopping
if args.stats:
    figure()
    suptitle( str(args.rundir) )
    data = genfromtxt( str(args.rundir)+"/stats.dat", dtype='float', comments='#', skip_footer=1 ).transpose() # transpose because we want the data column-wise
    maxdata = min( len(data[0])-1, args.tmax )
    # t	E/keV	V/keV	L/(m*kg*m/s)	P/(kg*m/s)
    subplot(321)
    xlabel("t/s")
    if len(data) > 8:
        ylabel(r"$\Delta\mathrm{E}/\mathrm{eV}$")
        plot( data[0][:maxdata], data[8][:maxdata]/N_Particles*1000 ) #factor 1000 needed, because stored in keV
    else:
        ylabel("E/eV")
        plot( data[0][:maxdata], data[1][:maxdata]/N_Particles*1000 ) #factor 1000 needed, because stored in keV

    subplot(322)
    xlabel("t/s")
    ylabel("V/eV")
    plot( data[0][:maxdata], data[2][:maxdata]/N_Particles*1000 )

    subplot(323)
    xlabel("t/s")
    ylabel("T/eV")
    if len(data) > 5:
        plot( data[0][:maxdata],  data[5][:maxdata]/N_Particles*1000 )
        plot( data[0][:maxdata], (data[1][:maxdata]-data[2][:maxdata])/N_Particles*1000 )
    else:
        plot( data[0][:maxdata], (data[1][:maxdata]-data[2][:maxdata])/N_Particles*1000 )

    subplot(324)
    xlabel("t/s")
    ylabel("P/(kg m/s)")
    plot( data[0][:maxdata], data[4][:maxdata]/N_Particles*1000 )

    # support for older stats.dat-files
    if len(data) > 7:
        subplot(325)
        xlabel("t/s")
        ylabel("Te/eV")
        plot( data[0][:maxdata], data[6][:maxdata]/N_Particles*1000 )

        subplot(326)
        xlabel("t/s")
        ylabel("Ti/keV")
        plot( data[0][:maxdata], data[7][:maxdata]/N_Particles*1000 )

# Last configuration in Simulation
if args.lastconfig:
    figure()
    if oldSimdataAvailable:
        data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )
        plot( data[-1], data[-2], "bo" ) # last ist x, because we skipped the last line, because this plot shall also work for non finished simulations
    else:
        dataEons = genfromtxt( str(args.rundir)+"/Electrons.dat", dtype='float', comments='#', skip_footer=1 )
        dataIons = genfromtxt( str(args.rundir)+"/Ions.dat",      dtype='float', comments='#', skip_footer=1 )
        plot( dataEons[-1], dataEons[-2], "bo" )
        plot( dataIons[-1], dataIons[-2], "ro" )

# Animation
if args.particlemovement:
    DPI   = matplotlib.rcParams['figure.dpi']

    newProgramOutputAvailable = False
    if os.path.isfile( str(args.rundir)+"/Electrons.dat" ):
        newProgramOutputAvailable = True
        dataEons = genfromtxt( str(args.rundir)+"/Electrons.dat", dtype='float', comments='#', skip_footer=1 )
        dataIons = genfromtxt( str(args.rundir)+"/Ions.dat",      dtype='float', comments='#', skip_footer=1 )
        if len(dataEons) > 0:
            if len(dataEons.shape) > 1:
                NEons = len(dataEons[0])
            else:
                NEons = 1
            dataEons = dataEons.reshape( (len(dataEons), NEons) ) # this only changes something if NEons == 1. Then it changes 1D to 2D array
        else:
            NEons = 0
        if len(dataIons) > 0:
            if len(dataIons.shape) > 1:
                NIons = NIons = len(dataIons[0])
            else:
                NIons = 1
            dataIons = dataIons.reshape( (len(dataIons), NIons) )
        else:
            NIons = 0
        print NIons,NEons
    if oldSimdataAvailable:
        data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )

    if oldSimdataAvailable and newProgramOutputAvailable:
        fig   = figure( figsize = ( 1200./DPI, 600./DPI ) )
        ax = subplot( 121, xlim=( 0, ceil(amax(data[0::2])) ),     ylim=( 0, ceil(amax(data[1::2])) ) )
        scatter, = ax.plot( data[0], data[1], "bo" )
        ax = subplot( 122, xlim=( 0, ceil(amax(dataEons[0::2])) ), ylim=( 0, ceil(amax(dataEons[1::2])) ) )
        if NEons > 0:
            scatterEons, = ax.plot( dataEons[0], dataEons[1], "bo" )
        if NIons > 0:
            scatterIons, = ax.plot( dataIons[0], dataIons[1], "ro" )
    elif oldSimdataAvailable:
        fig = figure( figsize = ( 800./DPI, 600./DPI ) )
        ax = subplot( 111, xlim=( 0, ceil(amax(data[0::2])) ),     ylim=( 0, ceil(amax(data[1::2])) ) )
        scatter, = ax.plot( data[0], data[1], "bo" )
    elif newProgramOutputAvailable:
        fig = figure( figsize = ( 800./DPI, 600./DPI ) )
        ax = subplot( 111, xlim=( 0, ceil(amax(dataEons[0::2])) ), ylim=( 0, ceil(amax(dataEons[1::2])) ) )
        if NEons > 0:
            scatterEons, = ax.plot( dataEons[0], dataEons[1], "bo" )
        if NIons > 0:
            scatterIons, = ax.plot( dataIons[0], dataIons[1], "ro" )

    subplots_adjust(bottom=0.1, left=.05, right=.95, top=.90, hspace=.35)
    suptitle( str(args.rundir) )

    # Animation control stuff
    clear = True
    paused = False
    currentFrame = 0
    def init():
        if oldSimdataAvailable:
            scatter.set_data([],[])
        if newProgramOutputAvailable:
            if NEons > 0:
                scatterEons.set_data([],[])
            if NIons > 0:
                scatterIons.set_data([],[])

        if oldSimdataAvailable and newProgramOutputAvailable:
            if NIons > 0 and NEons > 0:
                return (scatter,) + (scatterEons,) + (scatterIons,)
            elif NIons > 0:
                return (scatter,) + (scatterIons,)
            else:
                return (scatter,) + (scatterEons,)
        elif newProgramOutputAvailable:
            if NIons > 0 and NEons > 0:
                return (scatterEons,) + (scatterIons,)
            elif NIons > 0:
                return (scatterIons,)
            else:
                return (scatterEons,)
        else:
            return scatter,
    def animate(i):
        global currentFrame
        if not paused:
            currentFrame = (currentFrame + 1) % numberOfFrames
        if clear:
            if oldSimdataAvailable:
                scatter.set_data( data[2*currentFrame], data[2*currentFrame+1] )
            if newProgramOutputAvailable:
                if NEons > 0:
                    scatterEons.set_data( dataEons[2*currentFrame], dataEons[2*currentFrame+1] )
                if NIons > 0:
                    scatterIons.set_data( dataIons[2*currentFrame], dataIons[2*currentFrame+1] )
        else: # if not clear, then draw all particles from t=0 to current t which will form a kind of density plot
            if oldSimdataAvailable:
                a = data[0]
                b = data[1]
            if newProgramOutputAvailable:
                if NEons > 0:
                    ae = dataEons[0]
                    be = dataEons[1]
                if NIons > 0: 
                    ai = dataIons[0]
                    bi = dataIons[1]
            for k in range(1,currentFrame):
                if oldSimdataAvailable:
                    a = np.concatenate( (a, data[2*k]  ) )
                    b = np.concatenate( (b, data[2*k+1]) )
                if newProgramOutputAvailable:
                    if NEons > 0:
                        ae = np.concatenate( (ae, dataEons[2*k]  ) )
                        be = np.concatenate( (be, dataEons[2*k+1]) )
                    if NIons > 0:
                        ai = np.concatenate( (ai, dataIons[2*k]  ) )
                        bi = np.concatenate( (bi, dataIons[2*k+1]) )
            if oldSimdataAvailable:
                scatter.set_data( a,b )
            if newProgramOutputAvailable:
                if NEons > 0:
                    scatterEons.set_data( ae,be )
                if NIons > 0:
                    scatterIons.set_data( ai,bi )

        if oldSimdataAvailable and newProgramOutputAvailable:
            if NIons > 0 and NEons > 0:
                return (scatter,) + (scatterEons,) + (scatterIons,)
            elif NIons > 0:
                return (scatter,) + (scatterIons,)
            else:
                return (scatter,) + (scatterEons,)
        elif newProgramOutputAvailable:
            if NIons > 0 and NEons > 0:
                return (scatterEons,) + (scatterIons,)
            elif NIons > 0:
                return (scatterIons,)
            else:
                return (scatterEons,)
        else:
            return scatter,
    # call the animator.  blit=True means only re-draw the parts that have changed.

    if newProgramOutputAvailable:
        numberOfFrames = dataEons.shape[0] / 2
    elif oldSimdataAvailable:
        numberOfFrames = data.shape[0] / 2

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=numberOfFrames, interval=20, blit=True)
    if args.saveani:
        anim.save( str(args.rundir)+"/Simulation2D60fps.mp4", fps=60, dpi=(DPI/800.*1920.), extra_args=['-vcodec', 'libx264'])

    def key_analyzer(event):
        global clear, paused, currentFrame
        if (event.key == 'c'):
            clear = not clear
        if (event.key == ' ' or event.key == 'p'):
            paused = not paused
        if (event.key == '+'):
            currentFrame += 1
        if (event.key == '-'):
            currentFrame -= 1


    connect('key_press_event', key_analyzer)

if args.density:
    if oldSimdataAvailable:
        data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )
    else:
        dataEons = genfromtxt( str(args.rundir)+"/Electrons.dat", dtype='float', comments='#', skip_footer=1 )
        dataIons = genfromtxt( str(args.rundir)+"/Ions.dat",      dtype='float', comments='#', skip_footer=1 )
        if len(dataEons) > 0 and len(dataIons) > 0:
            data     = concatenate( (dataEons, dataIons) )
        elif len(dataEons) > 0:
            data = dataEons
        else:
            data = dataIons

    fig   = figure()
    DPI   = fig.get_dpi()
    # Estimate the 2D histogram
    nbins = max(32, int(sqrt(len(data)/500.)))
    print "N bins: ", nbins # so we have N data then we have sqrt(N) bins, which can all hold sqrt(N) datapoints if it was uniformly distributed, therefore making linear thigns visible
    y = data[1::2].reshape(-1)
    x = data[0::2].reshape(-1)
    if len(y) > len(x):
        print "shouldn't happen!"
    elif len(y) < len(x):
        x = data[0:-2:2].reshape(-1)
    maxdata = min( len(x), args.tmax )
    x = x[:maxdata]
    y = y[:maxdata]
    xmax = ceil(x.max()+0.1)
    ymax = ceil(y.max()+0.1)

    H, xedges, yedges = histogram2d( x,y, bins=nbins, range=[[0, xmax], [0, ymax]] )
    H = flipud(rot90(H)) # H needs to be rotated and flipped
    print "Hsum: ", H.sum(), " H.min:", H.min()," H.max:", H.max()
    H = H / H.sum()

    # Plot 2D histogram using pcolor
    axis( [ 0,xmax, 0,ymax ] )
    pcolormesh( xedges, yedges, H )
    cbar = colorbar()
    #cbar.ax.set_ylabel('Teilchenzahl')

    #tight_layout() #doesn't work with colobar
    fig.savefig( args.rundir+'_density.pdf', format='PDF')

show()
exit()
