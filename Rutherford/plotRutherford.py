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
parser.add_argument("-s", "--stats", help="", action="store_true")
parser.add_argument("-a", "--particlemovement", help="", action="store_true")
parser.add_argument("--tmax", help="", type=float)
args = parser.parse_args()

# loadtext/genfromtxt always loads the file row-wise in an array data[0] being the first row

if args.stats:
    figure()
    data = genfromtxt( str(args.rundir)+"/traj.dat", dtype='float', comments='#', skip_footer=1 ).transpose() # transpose because we want the data column-wise
    maxdata = min( len(data[0])-1, args.tmax )
    #t	rx	ry	rz	px	py	pz	r	theta	(Et-E)/E	(Lt-L)/L	theta analyt.	(th-th_anal)/th_anal
    #0  1   2   3   4   5   6   7     8         9          10            11                  12
    subplot(221)
    xlabel("x/m")
    ylabel("y/m")
    title("Trajectory")
    plot( data[1][:maxdata], data[2][:maxdata] )

    subplot(222)
    xlabel("t/s")
    ylabel(r"$\Delta E / E$")
    plot( data[0][:maxdata], data[9][:maxdata] )

    subplot(223)
    xlabel("t/s")
    ylabel(r"$\Delta L / L$")
    plot( data[0][:maxdata], data[10][:maxdata] )

    subplot(224)
    xlabel("t/s")
    ylabel(r"$\Delta \theta / \theta$")
    plot( data[0][:maxdata], data[12][:maxdata] )

show()
exit()
