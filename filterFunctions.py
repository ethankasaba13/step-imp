# Description: Contains functions used to design filters

# imports
import numpy as np
from pyaedt import Hfss

# Functions for finding g values
def findgs(order, ripple, filtertype):
    if filtertype == 'butter':
        gs = findbuttergs(order)
    elif filtertype == 'cheby1':
        gs = findcheby1gs(order, ripple)
    else:
        print('Error: filter type ' + filtertype + ' is not supported')

        gs = []
    return gs

def findcheby1gs(order, ripple):
    beta = np.log(1/np.tanh(ripple / 17.37)) # phase constant/beta value for Chebyshev response, Hong p. 41
    gamma = np.sinh(beta / (2 * order)) # propagation constant/gamma value for Chebyshev response, Hong p. 41
    gs = np.zeros(order + 2) # g values for Chebyshev response
    for i in range(order + 2):
        if i == 0:
            gs[i] = 1 # g0 value for Chebyshev response, see variables that should not be changed above
        elif i == 1:
            gs[i] = (2 / gamma) * np.sin(np.pi / (2 * order)) # g1 value for Chebyshev response, Hong p. 41
        elif i <= order:
            gs[i] = (1 / gs[i - 1]) * (4 * np.sin((2 * i - 1) * np.pi / (2 * order)) * np.sin((2 * i - 3) * np.pi / (2 * order))) / (gamma ** 2 + (np.sin(((i - 1) * np.pi) / order) ** 2)) # g n value for Chebyshev response, Hong p. 41 (3.26)
        else:
            if order % 2 == 0:
                gs[i] = (1 / np.tanh(beta / 4)) ** 2 # g n+1 value for Chebyshev response, Hong p. 41
            else:
                gs[i] = 1 # g n+1 value for Chebyshev response, Hong p. 41
    return gs

def findbuttergs(order):
    gs = np.zeros(order + 2)
    for i in range(order + 2):
        if i == 0:
            gs[i] = 1
        elif i <= order:
            gs[i] = 2 * np.sin((2 * i - 1) * np.pi / (2 * order)) # g n value for Butterworth response, Hong p. 41 (3.27)
        else:
            gs[i] = 1
    return gs

# Functions for finding component values
def inductancefromg(z0, g0, f, g): # Gets inductance from g value, Hong p. 114 (5.1)
    return (z0/g0) * (1/(2*np.pi*f)) * g
def capacitancefromg(z0, g0, f, g): # Gets capacitance from g value, Hong p. 114 (5.1)
    return (g0/z0) * (1/(2*np.pi*f)) * g

# Functions used in finding dimensions
def findulessthan2(ep, z): # finds ratio of width to height when w/h <= 2, done for inductor and source impedence
    A = (z / 60) * ((ep + 1) / 2) ** 0.5 + ((ep - 1) / (ep + 1)) * (0.23 + 0.11 / ep) # Used in next declaration, Hong p. 79
    return 8 * np.exp(A) / (np.exp(2 * A) - 2) # Hong p. 79 (4.10)

def findumorethan2(ep, z): # finds ratio of width to height when w/h >= 2, done for capacitor
    B = (60 * np.pi ** 2) / (z * np.sqrt(ep)) # Used in next declaration, Hong p. 79
    return (2 / np.pi) * (B - 1 - np.log(2 * B - 1) + ((ep - 1) / (2 * ep)) * (np.log(B - 1) + 0.39 - 0.61 / ep)) # Hong p. 79 (4.11)

def geteffepsilon(ep, u): # finds effective permittivity, Hong p. 78 (4.4)
    a = 1 + 1/49 * np.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) + 1/18.7 * np.log(1 + (u / 18.1) ** 3)
    b = 0.564 * ((ep - 0.9) / (ep + 3)) ** 0.053
    return ((ep + 1) / 2) + ((ep - 1) / 2) * ((1 + 10 / u) ** (-a * b))

def findlambda(cutoff, effepsilon): # find guided wavelength, Hong p. 77 (4.6b)
    return 300 / ((cutoff) * np.sqrt(effepsilon)) # guided wavelength (mm)

def getNewLengthL(lamd, cutoff, l, z0l): # finds length, Hong p. 114 (5.2) 
    return (lamd / (2 * np.pi)) * np.arcsin((cutoff * l) / z0l) 

def getNewLengthC(lamd, cutoff, c, z0c): # finds length, Hong p. 114 (5.2)
    return (lamd / (2 * np.pi)) * np.arcsin(cutoff * c * z0c)

# Functions for HFSS
def generateLines(hfss, lengths, widths, subheight, subwidth): # Generates high- and low-impedance lines in HFSS
    order = len(lengths) # order of filter
    lengthTracker = 0 # keeps track of length of filter
    if order < 3 or order > 9: # order must be between 3 and 9 inclusive
        print("Error: Order must be between 3 and 9 inclusive")
        return

    L1 = hfss.modeler.create_rectangle('XY', [0, subwidth/2 - widths[0]/2, subheight], [lengths[0], widths[0]], name='L1') # creates first low-impedance line
    L1.color = ('255', '0', '0')
    lengthTracker += lengths[0]

    C2 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[1]/2, subheight], [lengths[1], widths[1]], name='C2') # creates first high-impedance line
    C2.color = ('255', '0', '0')
    lengthTracker += lengths[1]

    L3 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[0]/2, subheight], [lengths[2], widths[0]], name='L3') # creates second low-impedance line
    L3.color = ('255', '0', '0')
    lengthTracker += lengths[2]

    if order >= 4:
        C4 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[1]/2, subheight], [lengths[3], widths[1]], name='C4') # creates second high-impedance line
        C4.color = ('255', '0', '0')
        lengthTracker += lengths[3]

    if order >= 5:
        L5 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[0]/2, subheight], [lengths[4], widths[0]], name='L5') # creates third low-impedance line
        L5.color = ('255', '0', '0')
        lengthTracker += lengths[4]

    if order >= 6:
        C6 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[1]/2, subheight], [lengths[5], widths[1]], name='C6') # creates third high-impedance line
        C6.color = ('255', '0', '0')
        lengthTracker += lengths[5]
    
    if order >= 7:
        L7 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[0]/2, subheight], [lengths[6], widths[0]], name='L7') # creates fourth low-impedance line
        L7.color = ('255', '0', '0')
        lengthTracker += lengths[6]

    if order >= 8:
        C8 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[1]/2, subheight], [lengths[7], widths[1]], name='C8') # creates fourth high-impedance line
        C8.color = ('255', '0', '0')
        lengthTracker += lengths[7]
    
    if order == 9:
        L9 = hfss.modeler.create_rectangle('XY', [lengthTracker, subwidth/2 - widths[0]/2, subheight], [lengths[8], widths[0]], name='L9') # creates fifth low-impedance line
        L9.color = ('255', '0', '0')