from numpy import *
from matplotlib.pyplot import * #xlabel,ylabel,plot,legend,show,step,connect,title,setp
from matplotlib import animation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tmax", help="", type=float)
args = parser.parse_args()

folders = [ "N25_Point_Periodic",
            "N25_Point_Periodic_Location_And_Force",
            "N25_Point_Reflecting",
            "N25_Point_wo_stopping_1e6" ]
folders = [ "N25_CIC_Periodic_dt1.5e-18",
            "N25_CIC_Periodic_Location_And_Force_dt1.5e-18",
            "N25_CIC_Reflecting_dt1.5e-18" ]
folders = [ "N25_CIC_Periodic_dt1.5e-18_wrong",
            "N25_CIC_Periodic_Location_And_Force_dt1.5e-18_wrong",
            "N25_CIC_Reflecting_dt1.5e-18_wrong" ]
folders = [ "N25_CIC_Periodic",
            "N25_CIC_Periodic_Location_And_Force",
            "N25_CIC_Reflecting" ]

folders = [ "N25_CIC_Periodic_dt1e-21",
            "N25_CIC_Periodic_dt1e-20",
            "N25_CIC_Periodic_dt1e-19",
            "N25_CIC_Periodic_dt1e-18" ]
pdfname = "N25_CIC_Periodic"

folders = [ "N25_CIC_Reflecting_dt1e-21",
            "N25_CIC_Reflecting_dt1e-20",
            "N25_CIC_Reflecting_dt1e-19",
            "N25_CIC_Reflecting_dt1e-18" ]
pdfname = "N25_CIC_Reflecting"

labels  = [r"$\Delta_t=10^{-21}\mathrm{s}$",
           r"$\Delta_t=10^{-20}\mathrm{s}$",
           r"$\Delta_t=10^{-19}\mathrm{s}$",
           r"$\Delta_t=10^{-18}\mathrm{s}$" ]


folders = [ "N25_CIC_Periodic_Location_And_Force_dt1e-19_r0.5cellWidth",
            "N25_CIC_Periodic_Location_And_Force_dt1e-19_r1cellWidth",
            "N25_CIC_Periodic_Location_And_Force_dt1e-19_r2cellWidth",
            "N25_CIC_Periodic_Location_And_Force_dt1e-19_r4cellWidth" ]
pdfname = "N25_CIC_Periodic_Force"

labels  = [r"$\mathrm{r}_\mathrm{CutOff}=0.5\cdot \mathrm{r}_\mathrm{Cloud}$",
           r"$\mathrm{r}_\mathrm{CutOff}=1.0\cdot \mathrm{r}_\mathrm{Cloud}$",
           r"$\mathrm{r}_\mathrm{CutOff}=2.0\cdot \mathrm{r}_\mathrm{Cloud}$",
           r"$\mathrm{r}_\mathrm{CutOff}=4.0\cdot \mathrm{r}_\mathrm{Cloud}$" ]


colors = ['r','g','c','k']
N_Particles = 25.
labelfontsize = 11

# loadtext always loads the file row-wise in an array data[0] being the first row
# E,T,L,V,P for random initial positions going into quilibrium through stopping
fig = figure( figsize=(14,6) )
for i in range(len(folders)):
    data = genfromtxt( folders[i]+"/stats.dat", dtype='float', comments='#', skip_footer=1).transpose()
    maxdata = min( len(data[0]), args.tmax )
    # transpose because we want the data column-wise
    subplot(131)
    xlabel("t/s")
    ylabel("<E>/eV")
    # t	E/keV	V/keV	L/(m*kg*m/s)	P/(kg*m/s)
    plot( data[0][:maxdata],data[1][:maxdata]/N_Particles*1000, label=labels[i], color=colors[i] )
    #*1000 because data unit is keV
    #legend( loc = 'center right', prop = {'size':labelfontsize} )

    subplot(132)
    xlabel("t/s")
    ylabel("<V>/eV")
    plot( data[0][:maxdata],data[2][:maxdata]/N_Particles*1000, label=labels[i], color=colors[i] )
    legend( loc = 'upper right', prop = {'size':labelfontsize} )
    #We can use only one legend, because we are sure, that the colors are given in the same order in all 3 plots

    subplot(133)
    xlabel("t/s")
    ylabel("<T>/eV")
    plot( data[0][:maxdata],data[5][:maxdata]/N_Particles*1000, label=labels[i], color=colors[i] )
    #legend( loc = 'upper left', prop = {'size':labelfontsize} )

    '''subplot(324)
    xlabel("t/s")
    ylabel("L/keV")
    plot( data[0],data[3], label=dir )
    legend( loc = 'upper left', prop = {'size':9} )

    subplot(325)
    xlabel("t/s")
    ylabel("P/keV")
    plot( data[0],data[4], label=dir )
    legend( loc = 'upper left', prop = {'size':9} )'''

tight_layout()
fig.savefig( pdfname+'.pdf', format='PDF')
show()
exit()


