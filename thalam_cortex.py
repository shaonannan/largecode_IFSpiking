import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio 
import os
#from connector_functions import gabor_probability
#from kernel_functions import gabor_kernel

# 12 * 10 * 6 * 10

n_pick = 8
g = 3.0

w = 0.8
phi = 0
gamma = 1  # Aspect ratio
sigma = 1
theta = 0

# Space parameters
dx = 0.3
lx = 3.0
dy = 0.3
ly = 3.0

xc = 0
yc = 0

x_values = np.arange(-lx/2, lx/2, dx)
y_values = np.arange(-ly/2, ly/2, dy)
# numbers of neurons
xnc = 20*24
ync = 6*24
nmax = xnc*ync
nlgnmax = 16


if os.path.exists("dlgnon.txt"):  
    os.remove("dlgnon.txt")  
if os.path.exists("dlgnoff.txt"):  
    os.remove("dlgnoff.txt")  

ft    = open("dlgnon.txt", 'w')  
ftoff = open("dlgnoff.txt", 'w')
align_dlgnon  = np.zeros((nmax,nlgnmax))
align_dlgnoff = np.zeros((nmax,nlgnmax))
def savetxt(filename,x):
    np.savetxt(filename,x,fmt=['%d ']*np.size(x),newline='\n')
def gabor_probability(x, y, sigma, gamma, phi, w, theta, xc=0, yc=0):

    """
    calculate the gabor function of x and y

    Returns value of the 2D Gabor function at x, y

    sigma: Controls the decay of the exponential term
    gamma: x:y proportionality factor, elongates the pattern
    phi: Phase of the overall pattern
    w: Frequency of the pattern
    theta: Rotates the whole pattern by the angle theta
    xc, yc : Linear translation
    """

    transforms_to_radians = np.pi / 180
    theta *= transforms_to_radians
    phi *= transforms_to_radians  # Transforms to radians

    # Translate
    x = x - xc
    y = y - yc

    # Rotate
    aux1 = np.cos(theta) * x + np.sin(theta) * y
    y = -np.sin(theta) * x + np.cos(theta) * y

    x = aux1

    # Function
    r = x**2 + (gamma * y) ** 2
    exp_part = np.exp(- r / (2 * sigma**2))
    cos_part = np.cos(2 * np.pi * w * x + phi)

    return exp_part * cos_part
def gabor_kernel(lx, dx, ly, dy, sigma, gamma, phi, w, theta, xc=0, yc=0):
    """
    Produces a gabor pattern. That is, the product of an exponential
    term and a sinusoidal term 
    
    Parameters
    -----------
    lx, ly : Wide of the spatial kernel in x and y respectively 
    dx, dy : Resolution in x and y respectively 
    sigma: Controls the decay of the exponential term
    gamma: x:y proportinality factor, elongates the pattern
    phi: Phase of the overall pattern 
    w: Frequency of the pattern 
    theta: Rotates the whole pattern by the angle theta
    
    """

    transforms_to_radians = np.pi / 180
    theta *= transforms_to_radians  # Transforms to radians
    phi *= transforms_to_radians

    x = np.arange(-lx/2, lx/2, dx)
    y = np.arange(-ly/2, ly/2, dy)

    # Translate
    x -= xc
    y -= yc

    X, Y = np.meshgrid(x, y)

    aux = np.cos(theta) * X + np.sin(theta) * Y
    Y = -np.sin(theta) * X + np.cos(theta) * Y

    X = aux

    exp_part = np.exp(-(X**2 + (gamma * Y)**2)/(2 * sigma**2))
    cos_part = np.cos(2 * np.pi * w * X + phi)

    return exp_part * cos_part


def GaborSample(x_value,y_value,sigma,gamma,phi,w,theta,n_pick,polarity):
    Z = np.zeros((x_values.size, y_values.size))
    L = np.zeros(Z.shape)
    xc = 0 
    yc = 0
    for x_index, x in enumerate(x_values):
        for y_index, y in enumerate(y_values):
            probability = polarity * gabor_probability(x, y, sigma, gamma, phi, w, theta, xc, yc)
            #print('proba', probability)
            counts = np.random.rand(n_pick) < probability
            #print('counts', counts)
            aux = np.sum(counts)  # Samples
            synaptic_weight = (g / n_pick) * aux
            L[x_index, y_index] = aux
            Z[x_index, y_index] = synaptic_weight
    return (L,Z)

sizex = x_values.size
sizey = y_values.size
countx = 0
county = 0
RFon  = {}
RFoff = {}
Bon   = {}
Boff  = {}
dlgnon = {}
dlgnoff = {}
# load theta and phase
theta  = np.loadtxt('theta.txt')
theta  = theta /np.pi * 180.0
phase  = np.loadtxt('phase.txt')
phase  = phase / 4.0 * 360.0
            
for yv in range(0,6*24,1):
    countx = 0
    for xv in range(0,20*24,1):
        # for lgnon/lgnoff index
        lgnon  = np.zeros((20*20,6*20))
        lgnoff = np.zeros((20*20,6*20)) 
        xc = np.round(xv/24.0*20.0)
        yc = np.round(yv/24.0*20.0)
        idxy   = xv + yv *20*24
        phi    = phase[idxy]
        the    = theta[idxy]
        #print(idxy,the)
        polarity = 1 # ON
        Lon,Zon   = GaborSample(x_values,y_values,sigma,gamma,phi,w,the,n_pick,polarity)
        polarity  = -1 # OFF
        Loff,Zoff = GaborSample(x_values,y_values,sigma,gamma,phi,w,the,n_pick,polarity)
        
        xs = xc - sizex/2
        xe = xc + sizex/2
        ys = yc - sizey/2
        ye = yc + sizey/2
        
        xstart = min(max(xs,0),20*20-1)
        xend   = min(max(xe,1),20*20)
        xids   = xstart - xs
        xrids  = int(xids)
        xide   = xe - xend
        xride  = int(sizex - xide)
        
        ystart = min(max(ys,0),6*20-1)
        yend   = min(max(ye,1),6*20)
        yids   = ystart - ys
        yrids  = int(yids)
        yide   = ye - yend
        yride  = int(sizey - yide)
        
        
        (xstart,xend,ystart,yend) = (int(xstart),int(xend),int(ystart),int(yend))
        RFon[idxy]   = Zon[xrids:xride,yrids:yride] 
        RFoff[idxy]  = Zoff[xrids:xride,yrids:yride] 
        
        Bon[idxy]   = Lon[xrids:xride,yrids:yride] 
        Boff[idxy]  = Loff[xrids:xride,yrids:yride] 
        lgnon[xstart:xend,ystart:yend] = Bon[idxy]
        x_y         = np.where(lgnon)
        dlgnon[idxy]  = x_y[:][0]+x_y[:][1]*20*20
        lgnoff[xstart:xend,ystart:yend]= Boff[idxy]
        x_y         = np.where(lgnoff)
        dlgnoff[idxy]  = x_y[:][0]+x_y[:][1]*20*20
        ttcount  = (county*3 + countx +1)
        
        #ft.write(str(np.reshape(dlgnon[idxy],(-1)))+'\n')
        for ion in range(np.size(dlgnon[idxy])):
            ft.write(str(dlgnon[idxy][ion])+' ')
            if ion<nlgnmax:
                align_dlgnon[idxy,ion] = dlgnon[idxy][ion]
        ft.write('\n')
        for ioff in range(np.size(dlgnoff[idxy])):
            ftoff.write(str(dlgnoff[idxy][ioff])+' ')
            if ioff<nlgnmax:
                align_dlgnoff[idxy,ioff] = dlgnoff[idxy][ioff]
        ftoff.write('\n')

        #plt.figure(ttcount)
        #plt.imshow(lgnon,cmap='jet')
#        plt.figure(2)
#        plt.subplot(6,12,ttcount)
#        plt.imshow(RFoff)
        countx +=1
        #print(countx,county,ttcount)
    county +=1
ft.close()      
ftoff.close()   
np.save("align_dlgnon.npy",align_dlgnon)   
np.save("align_dlgnoff.npy",align_dlgnoff)  

