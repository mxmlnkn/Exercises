from numpy import *
from matplotlib.pyplot import * #xlabel,ylabel,plot,legend,show,step,connect,title,setp
from matplotlib import animation
import argparse

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
folders = [ "N25_CIC_Reflecting_dt1.5e-18",
            "N25_CIC_Reflecting_dt1.5e-19",
            "N25_CIC_Reflecting_dt1.5e-20" ]
folders = [ "N25_CIC_Periodic_dt1.5e-18",
            "N25_CIC_Periodic_dt1.5e-19",
            "N25_CIC_Periodic_dt1.5e-20",
            "N25_CIC_Periodic_dt1.5e-21" ]

# loadtext always loads the file row-wise in an array data[0] being the first row
# E,T,L,V,P for random initial positions going into quilibrium through stopping
figure()
for dir in folders:
    data = genfromtxt( dir+"/stats.dat", dtype='float', comments='#', skip_footer=1).transpose()
    # transpose because we want the data column-wise
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
    plot( data[0],data[1]-data[2], label=dir )
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

show()
exit()


