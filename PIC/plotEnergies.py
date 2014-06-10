from numpy import *
from matplotlib.pyplot import * #xlabel,ylabel,plot,legend,show,step,connect,title,setp
from matplotlib import animation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("rundir", help="Directory which contains 'simOutput'")
parser.add_argument("-a", "--saveani", help="save animation of 2D Simulation", action="store_true")
args = parser.parse_args()

# loadtext/genfromtxt always loads the file row-wise in an array data[0] being the first row

# Energy Distributions for uniformly random start positions

figure()
data = genfromtxt( str(args.rundir)+"/CellEnergies.dat", dtype='float', comments='#' )
step( data[0],data[1], where='mid', label="Coulomb" )
step( data[0],data[2], where='mid', label="Spheres" )
step( data[0],data[3], where='mid', label="CIC" )
legend()

# E,T,L,V,P for random initial positions going into quilibrium through stopping

figure()
dir = str(args.rundir)
data = genfromtxt( str(args.rundir)+"/stats.dat", dtype='float', comments='#', skip_footer=1 ).transpose() # transpose because we want the data column-wise
# t	E/keV	V/keV	L/(m*kg*m/s)	P/(kg*m/s)
subplot(321)
xlabel("t/s")
ylabel("E/keV")
# t	E/keV	V/keV	L/(m*kg*m/s)	P/(kg*m/s)
plot( data[0],data[1], label=dir )
legend( loc = 'upper left', prop = {'size':9} )

subplot(322)
xlabel("t/s")
ylabel("V/keV")
plot( data[0],data[2], label=dir )
legend( loc = 'upper left', prop = {'size':9} )

subplot(323)
xlabel("t/s")
ylabel("T/keV")
plot( data[0],data[5], label=dir )
legend( loc = 'upper left', prop = {'size':9} )

subplot(324)
xlabel("t/s")
ylabel("L/keV")
plot( data[0],data[3], label=dir )
legend( loc = 'upper left', prop = {'size':9} )

subplot(325)
xlabel("t/s")
ylabel("P/keV")
plot( data[0],data[4], label=dir )
legend( loc = 'upper left', prop = {'size':9} )

# Last configuration in Simulation

data  = genfromtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#', skip_footer=1 )
figure()
plot( data[-2], data[-1], "bo" )

# Animation

fig   = figure()
DPI   = fig.get_dpi()
fig.set_size_inches( 800.0/float(DPI), 600.0/float(DPI) )
ax    = axes()#xlim=(0, 1), ylim=(0, 1))
scatter, = ax.plot( data[0], data[1], "bo" )
def init():
    scatter.set_data([],[])
    return scatter,
def animate(i):
    scatter.set_data(data[2*i],data[2*i+1])
    return scatter,
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(data[:,1])/2, interval=20, blit=True)
if args.saveani:
    anim.save( str(args.rundir)+"/Simulation2D60fps.mp4", fps=60, dpi=(DPI/800.*1920.), extra_args=['-vcodec', 'libx264'])

show()
exit()


