import netCDF4
import numpy as np
import matplotlib.pyplot as plt

import sys


# Open the NetCDF file
#filename = 'flora.nc'  # Replace with your NetCDF file path
filename = sys.argv[1]
nc = netCDF4.Dataset(filename, 'r')


var = nc.variables
dim = nc.dimensions

xro = var['XRO'][:]  # Read all data from 'XRO' variable
B = var['B'][:]  
try:
    dt = var['dt'][:]
except:
    # backward compatible
    dt = var['DT'][:]

mm = var['mm'][:]
flute3 = var['FLUTE3'][:]

Nt = dim['t'].size
Nz = dim['z'].size
Np = dim['psi'].size

xro2d = xro.reshape(Nt,Nz-2,Np-2)

B2d = B.reshape(Nz,Np)
jmax = np.argmax(B2d[:,0])
jmin = jmax + np.argmin(B2d[jmax:,0])
j0 = 0

dlogx = xro2d[-1]/xro2d[-2] - 1
gamma = np.mean(dlogx/dt) # MHD growth rate
print(f"Growth Rate {gamma} (m = {mm})")
print(f"Flute3 {flute3}")


