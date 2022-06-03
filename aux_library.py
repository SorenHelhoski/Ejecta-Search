import numpy as np
import matplotlib.pyplot as plt

#=============================================================

def Extract(file0):
    return np.array(list(eval(file0.readline().replace('\n',''))))

#=============================================================

def better_average(x):
    if len(x)==0:
        return None
    else:
        return np.average(x)

def better_std(x):
    if len(x)==0:
        return None
    else:
        return np.std(x)

def better_sum(x):
    if len(x)==0:
        return 0
    else:
        return sum(x)

def better_arctan(top,bottom):
    if top == 0 and bottom == 0 :
        return 0 # technically undefined
    if top == 0 and bottom > 0 :
        return 0
    if top == 0 and bottom < 0 :
        return np.pi
    if top > 0 and bottom == 0 :
        return np.pi/2
    if top < 0 and bottom == 0 :
        return 3*np.pi/2
    if bottom > 0 :
        return (np.arctan(top/bottom))%(2*np.pi)
    if bottom < 0 :
        return (np.arctan(top/bottom)+np.pi)%(2*np.pi)
    print('How did this happen?')

def SumIf(sum_list, cond_list, truths = 1):
    total = 0
    if type(truths) == 'int' or 'float':
        truths = [truths]
    for i in range(len(sum_list)):
        if cond_list[i] in truths:
            total += sum_list[i]
    return total

#=============================================================

def Basin_Size(d_imp,d_target,R_p,M_p,R_imp,v_imp):
    G=6.67e-11 #m,kg,s
    L=1.15
    R_p*=1000 ; R_imp*=1000 ; v_imp*=1000
    g = G*M_p/R_p**2
    return .5*1e-3*1.161*(d_imp/d_target)**(1.0/3.0)*(2.0*R_imp)**.78*v_imp**.44*g**-.22
    
#=============================================================

# Manual Histogram Class

'''
This class takes a x list and then bins it into a histogram,
and stores the non-zero bins and the respective frequencies in each

bins is the number of bins

bounds determine the upper and lower limit of counting frequencies
'''

class Bin:
    def __init__(self, x, y = [], weight=[], bins = 50, bounds = []):
        
        #set the y list if no value is given
        if y == []:
            for i in range(len(x)):
                y.append(1)
        #set the weights if no value is given
        if weight == []:
            for i in range(len(x)):
                weight.append(1)
        
        if x != []:
            x,y,weight = zip(*sorted(zip(x,y,weight)))
            if bounds == []:
                bounds = [min(x),max(x)]
        else:
            bounds = [0,1] # dummy values to keep from crashing when empty

        lower, upper = bounds
        bin_size = (upper - lower) / bins

        # inintialize certain parameters
        bottom = lower # lower bound of bin

        x_bin = [] # average of each
        x_err = [] # std of each bin

        y_bin = [] # average of second variable
        y_err = [] # std of each bin second variable

        f_bin = [] # frequency in each bin
        f_err = [] # error in the frequency

        x_ = [] # temporary list of each bin
        y_ = [] # temporary list of second variable
        f_ = [] # temporary list of weighted frequencies
    
        i = 0
        final = len(x)
        for each in range(len(x)):
            if x[each] < lower:
                i += 1
            if x[each] > upper:
                final = each
                break
        # loops over data
        while i < final:
            # append to temporary bin list if within the bounds
            if bottom <= x[i] <= bottom + bin_size:
                x_.append(x[i])
                y_.append(y[i])
                f_.append(weight[i])
                i += 1
            else:
                x_bin.append(better_average(x_))
                x_err.append(better_std(x_))
                y_bin.append(better_average(y_))
                y_err.append(better_std(y_))
                f_bin.append(better_sum(f_))
                f_err.append(better_sum(f_)**(1/2))
                x_ = [] ; y_ = [] ; f_ = []
                bottom += bin_size
        x_bin.append(better_average(x_))
        x_err.append(better_std(x_))
        y_bin.append(better_average(y_))
        y_err.append(better_std(y_))
        f_bin.append(better_sum(f_))
        f_err.append(better_sum(f_)**(1/2))
        
        c = 0 # total in each bin
        cmlt = [] # cumulative total in each bin
        for each in f_bin:
            cmlt.append(c)
            c += each

        c = 0 # total in each bin
        cmltr = [] # cumulative total in each bin from the right
        for i in range(len(f_bin),0,-1):
            cmltr.append(c)
            c += f_bin[i-1]
        cmltr.reverse()

        self.x_bin    = x_bin
        self.x_err    = x_err
        self.y_bin    = y_bin
        self.y_err    = y_err
        self.f_bin    = f_bin
        self.f_err    = f_err
        self.cmlt     = cmlt
        self.cmltr    = cmltr
        self.bin_size = bin_size

    def get_cmlt(self, direction = 'lr'): # returns the cumulative percent dist
        if direction == 'lr':
            cmlt_ = []
            _max = self.cmlt[-1]
            for each in self.cmlt:
                cmlt_.append(float(each)/_max)
            return cmlt_
        if direction == 'rl':
            cmlt_ = []
            _max = self.cmltr[0]
            for each in self.cmltr:
                cmlt_.append(float(each)/_max)
            return cmlt_
        else:
            print('direction must "lr" or "rl"')

    def get_x(self, factor=1): # returns the x and x_err
        dummy = []
        for each in self.x_bin:
            if each == None:
                pass
            else:
                dummy.append(float(each)*factor)
        dummy0 = []
        for each in self.x_err:
            if each == None:
                pass
            else:
                dummy0.append(float(each)*factor)
        return np.array(dummy), np.array(dummy0)

    def get_y(self, factor=1): # returns the y and y_err
        dummy = []
        for each in self.y_bin:
            if each == None:
                pass
            else:
                dummy.append(float(each)*factor)
        dummy0 = []
        for each in self.y_err:
            if each == None:
                pass
            else:
                dummy0.append(float(each)*factor)
        return np.array(dummy), np.array(dummy0)

    def get_f(self, factor=1): # returns the freq and f_err
        dummy = []
        for each in self.f_bin:
            if each == 0:
                pass
            else:
                dummy.append(float(each)*factor)
        dummy0 = []
        for each in self.f_err:
            if each == 0:
                pass
            else:
                dummy0.append(float(each)*factor)
        return np.array(dummy), np.array(dummy0)

#============================================================
#                 Sphereical Range Formula
#============================================================

def Coords(x,y,Dx,Dy,Dt,R,shape):
    velocity = (np.sqrt(Dx**2+Dy**2))/Dt
    if shape == 'DEFAULT':
        angle = better_arctan(Dy,Dx)
        r0 = R + y
        th0 = x/R
        return velocity, angle, r0, th0
    if shape == 'PLANET':
        angle = (better_arctan(Dy,Dx) + better_arctan(x,(y+R)))%(2*np.pi)
        r0 = np.sqrt(x**2 + (y+R)**2)
        th0 = better_arctan(x,(y+R))
        return velocity, angle, r0, th0
    

def Landing(initial_params, R, M, freq=10, L1 = np.inf):
    G = 6.67e-20 #km,kg,s
    velocity, angle, r0, th0 = initial_params
    t = [0] ; r = [r0] ; th = [th0]
    r_v = [velocity*np.sin(angle)] ; th_v = [velocity*np.cos(angle)/r0]

    if velocity > np.sqrt(2*G*M/r0):
        return 'escaped',None,None,None

    while abs(th[-1]-th[0]) < 2*np.pi :
        dt = (1/freq)*(r[-1]/R)**2
        t.append(t[-1] + dt)
        ra = r[-1]*th_v[-1]**2-G*M*r[-1]**-2
        tha = -2*th_v[-1]*r_v[-1]*r[-1]**-1
        r_v.append(r_v[-1]+ra*dt)
        th_v.append(th_v[-1]+tha*dt)
        r.append(r[-1]+r_v[-1]*dt)
        th.append(th[-1]+th_v[-1]*dt)
        if r[-1] < R:
            return t[-1], th[-1], R*th[-1],max(r)
        if r[-1] > L1:
            return 'distant',None,None,None # Goes up too high
    return 'orbit',None,None,None # Enters into an orbit







