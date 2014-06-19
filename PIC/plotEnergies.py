"""
Todo:
    - only average trajectory over certain time interval, which maybe can be given ...
"""

from numpy import *
from matplotlib.pyplot import * #xlabel,ylabel,plot,legend,show,step,connect,title,setp
from matplotlib import animation
import argparse

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

# loadtext/genfromtxt always loads the file row-wise in an array data[0] being the first row

# Energy Distributions for uniformly random start positions
if args.initenergies:
    figure()
    data = genfromtxt( str(args.rundir)+"/CellEnergies.dat", dtype='float', comments='#' )
    step( data[0],data[1], where='mid', label="Coulomb" )
    step( data[0],data[2], where='mid', label="Spheres" )
    step( data[0],data[3], where='mid', label="CIC" )
    legend()

# E,T,L,V,P for random initial positions going into quilibrium through stopping
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
    data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )
    figure()
    plot( data[-1], data[-2], "bo" ) # last ist x, because we skipped the last line, because this plot shall also work for non finished simulations

# Animation
if args.particlemovement:
    data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )
    fig   = figure()
    DPI   = fig.get_dpi()
    fig.set_size_inches( 800.0/float(DPI), 600.0/float(DPI) )
    ax    = axes()#xlim=(0, 1), ylim=(0, 1))
    scatter, = ax.plot( data[0], data[1], "bo" )
    title( str(args.rundir) )
    clear = True
    def init():
        scatter.set_data([],[])
        return scatter,
    def animate(i):
        if clear:
            scatter.set_data( data[2*i], data[2*i+1] )
        else:
            a = data[0]
            b = data[1]
            for k in range(1,i):
                a = np.concatenate( (a, data[2*k]  ) )
                b = np.concatenate( (b, data[2*k+1]) )
            scatter.set_data( a, b )
        return scatter,
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(data[:,1])/2, interval=20, blit=True)
    if args.saveani:
        anim.save( str(args.rundir)+"/Simulation2D60fps.mp4", fps=60, dpi=(DPI/800.*1920.), extra_args=['-vcodec', 'libx264'])

    def key_analyzer(event):
        global clear
        if (event.key == 'c'):
            print "clear was: ",clear
            clear = not clear

    connect('key_press_event', key_analyzer)

if args.density:
    data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )
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

    H, xedges, yedges = histogram2d( x,y, bins=nbins )
    H = flipud(rot90(H)) # H needs to be rotated and flipped
    print "Hsum: ", H.sum(), " H.min:", H.min()," H.max:", H.max()
    H = H / H.sum()

    # Plot 2D histogram using pcolor
    axis([0,int(x.max()+0.1), 0,int(y.max()+0.1)])
    pcolormesh(xedges,yedges,H)
    cbar = colorbar()
    #cbar.ax.set_ylabel('Teilchenzahl')

    #tight_layout() #doesn't work with colobar
    fig.savefig( args.rundir+'_density.pdf', format='PDF')

show()
exit()


