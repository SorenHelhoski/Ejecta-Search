import numpy as np
import matplotlib.pyplot as plt
import pySALEPlot as psp
from aux_library import better_arctan, Coords, Landing, SumIf
from math import floor
from tqdm import tqdm # you have to install this package

import os,glob
import sys

show_output = True
skips = 1

original = sys.stdout
supress = open(os.devnull,'w')
if show_output:
    supress = original



#============================================================
#                      Custom Parameters 
#============================================================

print('Extracting Ejecta.inp Parameters')

input_file = open('./Ejecta.inp','r')
N = None
param = {'NME':N,'FTS':N,'BN1':N,'BN2':N,'DN0':N,'DN1':N,'DN2':N,'DN3':N,'MPL':N,'RPL':N,'MTS':N,'LGI':N,'NCO':N,'ECO':N,'AST':N,'SIM':N}
for line in input_file:
    row = line.split(':')
    KEY = (row[0]).replace(' ','')
    VAL = (row[-1][0:-1]).replace(' ','')
    if KEY=='AST' or KEY=='SIM' or KEY=='NME':
        param[KEY] = str(VAL)
    elif KEY in param.keys():
        param[KEY] = eval(VAL)

name = param['NME'] # name of saved data file

time_f = param['FTS'] # Final time STEP

dens0 = param['DN0']*1e12;dens1 = param['DN1']*1e12;dens2 = param['DN2']*1e12;dens3 = param['DN3']*1e12 # densities of layers

M_p = param['MPL'] # Mass of Planet
R_p = param['RPL'] # Radius of Planet

mts = param['MTS'] # minimum timestep of balistic trajectory
lgi = param['LGI'] # max radius of range simulator

lat_long = [param['NCO'],param['ECO']]

ast_file = param['AST']
datafile = param['SIM']

dirname = '.'.format(name)

try:
    os.makedirs(dirname+'/Highlights_{}'.format(name))
except:
    pass
if show_output:
    for filename in glob.glob(dirname+'/Highlights_{}/*'.format(name)):
        os.remove(filename)

print('DONE\n')

#============================================================
#                Parameters from asteroid.inp
#============================================================

print('Extracting asteroid.inp Parameters')

# open input file
ast_input = open(ast_file,'r')
keywords = ['GRIDSPC', 'OBJRESH', 'OBJVEL', 'DTSAVE','LAYPOS','GRIDH']
ast_dict = {}
for line in ast_input:
    word  = line[0:16].replace(' ','')
    if word == 'S_TYPE':
        Type = line[54:-1].replace(' ','')
    value = '['+(line[54:-1].replace(' ','').replace(':',',')).replace('D','e')+']'
    if word in keywords:
        ast_dict[word] = eval(value)

spacing = ast_dict['GRIDSPC'][0] *.001     # (km)
R_imp = ast_dict['OBJRESH'][0] * spacing   # radius of impactor (km)
v_imp = ast_dict['OBJVEL'][0]  * -.001     # impact velcity (km/s)
dt   = ast_dict['DTSAVE'][0]               # save interval of sim (s)
max_x= ast_dict['GRIDH'][1] * spacing 

if Type == 'DEFAULT':
    BHigh = -1e16 ; BLow  = -1e16
    layers = ast_dict['LAYPOS']
    if len(layers) > 1:
        BHigh = (layers[-2]-layers[-1])*spacing
    if len(layers) > 2:
        BHigh = (layers[-3]-layers[-1])*spacing

if Type == 'PLANET':
    BHigh = 0 ; BLow  = 0
    layers = ast_dict['OBJRESH']
    if len(layers) > 2:
        BHigh = layers[2]*spacing
    if len(layers) > 3:
        BLow = layers[3]*spacing

dx   = spacing      # tracer spacing horizontal
dy   = spacing      # tracer spacing vertical

print('DONE\n')
   
#============================================================
#                 Pre-Ejecta Search Extraction
#============================================================
model=psp.opendatfile(datafile)
model.setScale('km')

#Find the peak pressures of all tracers
print('Extracting Peak Pressure')
step = model.readStep('TrP', model.nsteps-1)
pres = np.append(np.zeros(0), step.TrP/10**9)
print('DONE\n')


# extract the materials and weights
print('Extracting Material, Volume, and Mass')

try:
    step = model.readStep('Trd', time)
    dd = 0.001*np.append(np.zeros(0), step.Trd) # get the density in [kg/km^3]
    dens0,dens1,dens2,dens3 = 1,1,1,1           # nullify custom density
except:
    dd = np.ones(len(pres)) 


step = model.readStep('TrT', 0)
xx_0 = np.append(np.zeros(0), step.xmark)
yy_0 = np.append(np.zeros(0), step.ymark)

#get the materials:
d_imp = 1 ; d_target = 1
mtr = [] # materials of all tracers
vol = [] # volume of all tracers [km^3]
mas = [] # mass of all tracers [kg]
for i in range(len(yy_0)): # sort based on initial location
    if Type == 'PLANET':
        V = abs(2*np.pi*xx_0[i]*dx*dy)
    if Type == 'DEFAULT':
        V = abs(2*np.pi*xx_0[i]*dx*dy*((R_p+yy_0[i])/R_p)) # extra factor corrects for change in curvature
    vol.append(V)
    if yy_0[i] > 0:
        mtr.append(0)
        mas.append(V*dens0*dd[i])
        d_imp = dens0*dd[i]
    elif (yy_0[i] > BHigh and Type=='DEFAULT') or (((yy_0[i]+R_p)**2+xx_0[i]**2 > (BHigh)**2) and Type=='PLANET'):
        mtr.append(1)
        mas.append(V*dens1*dd[i])
        d_target = dens1*dd[i]
    elif (yy_0[i] > BLow and Type=='DEFAULT') or (((yy_0[i]+R_p)**2+xx_0[i]**2 > (BLow)**2) and Type=='PLANET'):
        mtr.append(2) 
        mas.append(V*dens2*dd[i])
    else:                # deepest layer
        mtr.append(3)
        mas.append(V*dens3*dd[i])

print('DONE\n')

# make sure Results.dat exists
open('{}/{}.dat'.format(dirname,name),'a').close()
filetest = open('{}/{}.dat'.format(dirname,name),'r')    
filetest.close() 
open('{}/{}.out'.format(dirname,name),'a').close()
filetest = open('{}/{}.out'.format(dirname,name),'r')    
filetest.close() 
        
#============================================================
#                          Search
#============================================================       

time_0 = int(floor(2*R_imp/(dt*v_imp)+1))

step = model.readStep('TrT', time_0-1)
xx_0 = np.append(np.zeros(0), step.xmark)
yy_0 = np.append(np.zeros(0), step.ymark)

# initialize lists
used=[]
material=[];volume=[];spread=[];mass=[];
launch_time=[];launch_polar=[];launch_length=[];launch_height=[];launch_velocity=[];launch_angle=[];
pressure=[];temperature=[];
landing_time=[];landing_polar=[];landing_length=[];max_altitude=[]

#initialize escape 

escaped = 0 ; x_e=[];y_e=[] ; V_e=0;M_e=0
distant = 0 ; x_d=[];y_d=[] ; V_d=0;M_d=0
orbit = 0   ; x_o=[];y_o=[] ; V_o=0;M_o=0
bound = 0   ; x_b=[];y_b=[] ; V_b=0;M_b=0

tracer_index_0 = range(0,len(xx_0),skips)
tracer_index = []

for each in tracer_index_0: 
    if yy_0[each] > -R_p/6:
        tracer_index.append(each)

#print(len(tracer_index))

print('============================')
print('Starting Ejecta Search...')
print('============================')

# Loop over Tracers for Considered Timesteps

if not show_output:
    Loop = tqdm(range(time_0 , time_f+1))
if show_output:
    Loop = range(time_0 , time_f+1)

for time in Loop:
    bound = len(x_b)
    escaped = len(x_e)
    distant = len(x_d)
    orbit = len(x_o)
    # Read the step
    try:
        sys.stdout = supress
        step = model.readStep('TrT', time)
        sys.stdout = original
    except:
        sys.stdout = original
        print('Timestep {} OUT OF RANGE'.format(time))
        time_f = time-1
        break

    temp = np.append(np.zeros(0), step.TrT)
    xx_1 = np.append(np.zeros(0), step.xmark)
    yy_1 = np.append(np.zeros(0), step.ymark)

    for each in tracer_index:
        # Compare with previous step
        if Type == 'DEFAULT':
            H = yy_1[each]
            if xx_1[each] > max_x:
                H = -1 # cancel the search if the distance is past the high resolution zone
        if Type == 'PLANET':
            H = np.sqrt(xx_1[each]**2+(yy_1[each]+R_p)**2)-R_p
        if H > 2*R_imp:
            D_x = xx_1[each] - xx_0[each] # Change in x
            D_y = yy_1[each] - yy_0[each] # Change in y

            velocity, angle, r0, th0 = Coords(xx_1[each],yy_1[each],D_x,D_y,dt,R_p,Type)
            initial_param = [velocity, angle, r0, th0]
            
            if not (0 < angle < np.pi/2) :
                # pass over tracers that are traveling down
                continue

            tracer_index.remove(each)

            landing = Landing(initial_param,R_p,M_p,freq=1/mts,L1=lgi)

            if landing[0] == 'escaped':
                x_e.append(xx_1[each]),y_e.append(yy_1[each])
                V_e+=vol[each] ; M_e+=mas[each]
                continue
            if landing[0] == 'distant':
                x_d.append(xx_1[each]),y_d.append(yy_1[each])
                V_d+=vol[each] ; M_d+=mas[each]
                continue
            if landing[0] == 'orbit':
                x_o.append(xx_1[each]),y_o.append(yy_1[each])
                V_o+=vol[each] ; M_o+=mas[each]
                continue
            
            x_b.append(xx_1[each]),y_b.append(yy_1[each])
            V_b+=vol[each] ; M_b+=mas[each]

            used.append(each) # indexes of previously used tracers
            launch_time.append(time*dt)
            launch_polar.append(th0)
            launch_length.append(R_p*th0)
            launch_height.append(r0-R_p)
            launch_velocity.append(velocity)
            launch_angle.append(angle)
            temperature.append(temp[each])

            landing_time.append(landing[0]+time*dt)
            landing_polar.append(landing[1])
            landing_length.append(landing[2])
            max_altitude.append(landing[3]-R_p)
    xx_0, yy_0 = xx_1, yy_1 # makes the current step the previous step

    #====================
    # Plot the Progress
    #====================
    if show_output:
        fig=plt.figure(figsize=(8,4))
        ax=fig.add_subplot(111,aspect='equal')
        ax.set_xlabel('r [km]')
        ax.set_ylabel('z [km]')
        if Type == 'DEFAULT':
            ax.set_xlim([0,20*R_imp])
            ax.set_ylim([-2*R_imp,6*R_imp])
        if Type == 'PLANET':
            ax.set_xlim([0,R_p/2])
            ax.set_ylim([-R_p/4,5*R_imp])
        p1=ax.pcolormesh(model.x,model.y,step.mat, cmap='Oranges',vmin=1,vmax=model.nmat+1)

        ax.set_title('{: 5.2f} s'.format(step.time))
        ax.scatter(x_b[bound:-1],y_b[bound:-1],c='black',s=3,linewidths=0)#,label='will land')
        ax.scatter(x_e[bound:-1],y_e[bound:-1],c='red',  s=3,linewidths=0)#,label='will escape')
        ax.scatter(x_d[bound:-1],y_d[bound:-1],c='blue' ,s=3,linewidths=0)#,label='will travel too far')
        ax.scatter(x_o[bound:-1],y_o[bound:-1],c='gold', s=3,linewidths=0)#,label='will go into orbit')
        #ax.legend()
        fig.savefig('{}/Highlights_{}/MatTmp-{:05d}.png'.format(dirname,name,time,dpi=300))
        plt.close()

print('COMPLETED SEARCH\n')

for each in used:
    pressure.append(pres[each])
    material.append(mtr[each])
    volume.append(vol[each])
    mass.append(mas[each])

vol_array = np.array(volume)
polar_array = np.array(landing_polar)
for each in range(len(volume)):
    spread.append(R_p*np.sin(abs(polar_array[each])))

bound = len(x_b)
escaped = len(x_e)
distant = len(x_d)
orbit = len(x_o)

output_string = '\n'

output_string += ('Name of Results : {} \n'.format(name))
output_string += ('jdata.dat    : {} \n'.format(datafile))
output_string += ('asteroid.inp : {} \n\n'.format(ast_file))

output_string += ('Final Ejecta Launch at {:.2e} seconds\n'.format(max(launch_time)))
output_string += ('            Timestep = {}\n\n'.format(int(max(launch_time)/dt)))

if int(max(launch_time)/dt) == time_f:
    output_string += ('WARNING: Ejecta curtain is still ejecting material \nConsider changing the Final Time Step of the search\n\n')

output_string += ('                                Tracers       Volume          Mass\n')
output_string += ('----------------------------------------------------------------------------\n')

output_string += ('Planet                        : ------- : {:.3e} km^3 : {:.3e} kg\n'.format(4*np.pi*R_p**3/3,M_p))
output_string += ('Simulation                    : {:07d} : {:.3e} km^3 : {:.3e} kg\n\n'.format(len(xx_1),sum(vol),sum(mas)))

output_string += ('Total Ejecta                  : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(bound+escaped+distant+orbit,V_b+V_e+V_d+V_o,M_b+M_e+V_d+M_o)) 

output_string += ('      Landed  (on surface)    : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(bound,V_b,M_b))
output_string += ('      Escaped (to infinity)   : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(escaped,V_e,M_e))
output_string += ('      Distant (out of bounds) : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(distant,V_d,M_d))
output_string += ('      In Orbit (elliptical)   : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(orbit,V_o,M_o))
output_string += ('\n')

output_string += ('Of the Landed Ejecta... \n') 

output_string += ('      Impactor                : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(material.count(0),SumIf(volume,material,truths=0),SumIf(mass,material,truths=0)))
output_string += ('      Top Layer               : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(material.count(1),SumIf(volume,material,truths=1),SumIf(mass,material,truths=1)))
output_string += ('      Middle Layer            : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(material.count(2),SumIf(volume,material,truths=2),SumIf(mass,material,truths=2)))
output_string += ('      Bottom Layer            : {:07d} : {:.3e} km^3 : {:.3e} kg\n'.format(material.count(3),SumIf(volume,material,truths=3),SumIf(mass,material,truths=3)))
output_string += ('\n')

print(output_string)
      
print('Writing Output File')
try:
    fileO = open('{}/{}.out'.format(dirname,name),'x')
except:
    fileO = open('{}/{}.out'.format(dirname,name),'r+')
    fileO.truncate(0)

fileO.write(output_string)

fileO.close()

try:
    file0 = open('{}/{}.dat'.format(dirname,name),'x')
except:
    file0 = open('{}/{}.dat'.format(dirname,name),'r+')
    file0.truncate(0)

prop = [time_0, time_f, dt, R_p, M_p, R_imp, v_imp, lat_long]

file0.write(str(prop)+'\n')

file0.write(str(used)+'\n')

file0.write(str(material)+'\n')
file0.write(str(volume)+'\n')
file0.write(str(spread)+'\n')
file0.write(str(mass)+'\n')

file0.write(str(launch_time)+'\n')
file0.write(str(launch_polar)+'\n')
file0.write(str(launch_length)+'\n')
file0.write(str(launch_height)+'\n')
file0.write(str(launch_velocity)+'\n')
file0.write(str(launch_angle)+'\n')

file0.write(str(pressure)+'\n')
file0.write(str(temperature)+'\n')

file0.write(str(landing_time)+'\n')
file0.write(str(landing_polar)+'\n')
file0.write(str(landing_length)+'\n')
file0.write(str(max_altitude)+'\n')

file0.close()
print('DONE\n')
      
