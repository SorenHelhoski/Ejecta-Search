import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from aux_library import Bin, Extract
from scipy.optimize import curve_fit as fit

input_file = open('./Ejecta.inp','r')
param = {'NME':None}
for line in input_file:
    row = line.split(':')
    KEY = (row[0]).replace(' ','')
    VAL = (row[-1][0:-1]).replace(' ','')
    if KEY=='NME':
        param[KEY] = str(VAL)

name = param['NME'] # name of saved data file


dirname = '/home/shelhosk/Desktop/Ejecta_Output/Bulk_{}'.format(name)
psp.mkdir_p(dirname)      

print('Opening data file...')
# open tracer file
file0    = open('/home/shelhosk/Desktop/Ejecta_Output/{}.dat'.format(name),'r')
time_0, time_f, save_step, R_p, M_p, R_imp, v_imp, coords = eval(file0.readline().replace('\n',''))
used = Extract(file0)
material = Extract(file0)
volume   = Extract(file0)
spread   = Extract(file0)
mass     = Extract(file0)
launch_time     = Extract(file0)
launch_polar    = Extract(file0)
launch_length   = Extract(file0)
launch_height   = Extract(file0)
launch_velocity = Extract(file0)
launch_angle    = Extract(file0)
pressure    = Extract(file0)
temperature = Extract(file0)
landing_time   = Extract(file0)
landing_polar  = Extract(file0)
landing_length = Extract(file0)
max_altitude   = Extract(file0)
print('DONE\n')

#------------------------------------------------------------
#                     Graphs
#------------------------------------------------------------
def Graph(nx,ny,ref=False):
    dirvar = '{}/'.format(dirname)+nx
    psp.mkdir_p(dirvar) 
    name = ny+' v. '+nx
    x=eval(nx.replace(' ','_').lower()) ; y=eval(ny.replace(' ','_').lower())
    fig = plt.figure(figsize=(12, 6)) ; ax=fig.add_subplot(111)
    ax.set_title(name+' : [km][kg][s][rad][GPa][K]')
    ax.set_xlabel(nx) ; ax.set_ylabel(ny)
    plt.scatter(x,y,s=1,linewidths=0.01)
    if ref:
        ax.plot([0,0],[min(y),max(y)],label='Impact')
        ax.plot([np.pi*R_p/2,np.pi*R_p/2],[min(y),max(y)],ls=':',label='Demipode')
        if max(x)>np.pi*R_p/2:
            ax.plot([np.pi*R_p,np.pi*R_p],[min(y),max(y)],ls='--',label='Antipode')
    ax.grid(True,ls='--',zorder=-15,alpha=.1); ax.legend()
    fig.savefig('{}/{}.png'.format(dirvar,nx+' v. '+ny))
    print('Saved: {}.png'.format(name))

def Hist(nx,nw,ref=False):
    dirvar = '{}/Hist/'.format(dirname)+nw
    psp.mkdir_p(dirvar)
    name = 'Histogram of '+nx
    x=eval(nx.replace(' ','_').lower())
    w=eval(nw.replace(' ','_').lower())
    fig = plt.figure(figsize=(12, 6)) ; ax=fig.add_subplot(111)
    ax.set_title(name+' : [km][kg][s][rad][GPa][K]')
    ax.set_xlabel(nx) ; ax.set_ylabel('Total {}'.format(nw))
    ax.hist(x,weights=w,bins=100)
    if ref:
        ax.plot([0,0],[0,sum(w)/(max(x)-min(x))],label='Impact')
        ax.plot([np.pi*R_p/2,np.pi*R_p/2],[0,sum(w)/(max(x)-min(x))],ls=':',label='Demipode')
        if max(x)>np.pi*R_p/2:
            ax.plot([np.pi*R_p,np.pi*R_p],[0,sum(w)],ls='--',label='Antipode')
    ax.grid(True,ls='--',zorder=-15,alpha=.1); ax.legend()
    fig.savefig('{}/{}.png'.format(dirvar,nx))
    print('Saved: {}.png'.format(name))


Graph('Launch Length','Launch Time')
Graph('Launch Length','Launch Polar')
Graph('Launch Length','Launch Height')
Graph('Launch Length','Launch Velocity')
Graph('Launch Length','Launch Angle')
Graph('Launch Length','Landing Length')

Graph('Launch Velocity','Launch Angle')
Graph('Launch Velocity','Landing Length')
Graph('Launch Velocity','Max Altitude')
Graph('Launch Velocity','Landing Time')

Graph('Landing Length','Landing Time',ref=True)
Graph('Landing Length','Landing Polar',ref=True)
Graph('Landing Length','Max Altitude',ref=True)

Graph('Launch Length','Max Altitude')
Graph('Launch Velocity','Max Altitude')
Graph('Launch Angle','Max Altitude')

Graph('Pressure','Temperature')
Graph('Pressure','Launch Velocity')
Graph('Pressure','Material')

Graph('Launch Time','Landing Time')
Graph('Launch Time','Launch Length')

Graph('Landing Length','Spread',ref=True)
Graph('Launch Length','Volume')
Graph('Launch Length','Mass')
Graph('Volume','Mass')

Graph('Material','Mass')
Graph('Material','Volume')
Graph('Material','Landing Length')

Hist('Launch Time','Mass')
Hist('Launch Length','Mass')
Hist('Launch Height','Mass')
Hist('Landing Time','Mass')
Hist('Landing Length','Mass',ref=True)

Hist('Launch Time','Volume')
Hist('Launch Length','Volume')
Hist('Launch Height','Volume')
Hist('Landing Time','Volume')
Hist('Landing Length','Volume',ref=True)

Hist('Max Altitude','Mass')
Hist('Max Altitude','Volume')

Hist('Pressure','Mass')
Hist('Pressure','Volume')
Hist('Temperature','Mass')
Hist('Temperature','Volume')





