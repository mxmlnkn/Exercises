from numpy import *
from matplotlib.pyplot import * #xlabel,ylabel,plot,legend,show,step,connect,title,setp
from matplotlib import animation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("rundir", help="Directory which contains 'simOutput'")
args = parser.parse_args()

# loadtext always loads the file row-wise in an array data[0] being the first row

# Energy Distributions for unformly random start positions

figure()
data = loadtxt( str(args.rundir)+"/CellEnergies.dat", dtype='float', comments='#')
step( data[0],data[1], where='mid', label="Coulomb" )
step( data[0],data[2], where='mid', label="Spheres" )
step( data[0],data[3], where='mid', label="CIC" )
legend()

# E,T,L,V,P for random initial positions going into quilibrium through stopping

figure()
data = loadtxt( str(args.rundir)+"/stats.dat", dtype='float', comments='#').transpose() # transpose because we want the data column-wise
# t	E/keV	V/keV	L/(m*kg*m/s)	P/(kg*m/s)
plot( data[0],data[1]/data[1][-1], label="E/%e" % (data[1][-1]) )
plot( data[0],data[2]/data[2][-1], label="V/%e" % (data[2][-1]) )
plot( data[0],data[3]/data[3][-1], label="L/%e" % (data[3][-1]) )
plot( data[0],data[4]/data[4][-1], label="P/%e" % (data[4][-1]) )
legend()

# Last configuration in Simulation

data  = loadtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#' )
figure()
plot( data[-2], data[-1], "bo" )

# Animation

data  = loadtxt( str(args.rundir)+"/simdata.dat", dtype='float', comments='#' )
fig   = figure()
ax    = axes(xlim=(0, 1), ylim=(0, 1))
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

show()
exit()


