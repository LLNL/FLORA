import netCDF4
import numpy as np
import matplotlib.pyplot as plt

import sys

class FloraReader:

    def __init__(self,filename):
        self.nc = netCDF4.Dataset(filename, 'r')
        self.filename = filename
        self.loadData()

    def loadData(self):
        
        nc = self.nc
        var = nc.variables
        dim = nc.dimensions
        
        self.xro = var['XRO'][:]  # Read all data from 'XRO' variable
        self.B = var['B'][:]  
        self.R = var['R'][:]  
        self.bvac = var['BVAC'][:]  
        self.rho = var['RHO'][:]
        self.pperp = var['PPERP'][:]
        self.ppar = var['PPAR'][:]
        self.pperpe = var['PPERPE'][:]
        
        self.dt = var['dt'][:]
        self.mm = var['mm'][:]
        self.flute3 = var['FLUTE3'][:]
        self.zmax = var['zmax'][:]
        
        self.tenergy = var['tenergy'][:]  # Read all data from 'XRO' variable
        self.enkin = var['enkin'][:]  # Read all data from 'XRO' variable
        self.enpot = var['enpot'][:]  # Read all data from 'XRO' variable
        Nt = dim['t'].size
        Nz = dim['z'].size
        Np = dim['psi'].size
        
        self.xro2d = self.xro.reshape(Nt,Nz-2,Np-2)
        
        B2d = self.B.reshape(Nz,Np)
        jmax = np.argmax(B2d[:,0])
        jmin = jmax + np.argmin(B2d[jmax:,0])
        j0 = 0
        
        self.z_axis = np.linspace(0,self.zmax,Nz) # later can add coordinate stretching
        self.psi_axis = np.linspace(0,1,Np)
        self.zax, self.pax = np.meshgrid(self.z_axis,self.psi_axis)
        
#        dlogx = xro2d[-1]/xro2d[-2] - 1
#        gamma = np.mean(dlogx/self.dt) # MHD growth rate


    def plotFlute3(self):

        plt.figure()
        plt.plot(self.psi_axis, self.flute3)
        plt.title(self.filename)
        plt.grid()

    def plotXRO1D(self, t = -1):

        fig,axs = plt.subplots(2,2, figsize=(12,5))
        zax = self.z_axis[1:-1]
        pax = self.psi_axis[1:-1]
        xi = self.xro2d[t]

        B = self.B[1:-1, 1:-1]
        R = self.R[1:-1, 1:-1]
        rBxi = R*B*xi

        axs[0,0].plot(pax, xi[0])
        axs[1,0].plot(pax, rBxi[0])
        axs[0,1].plot(zax, xi[:,0])
        axs[1,1].plot(zax, rBxi[:,0])

        axs[0,0].set_ylabel(r"$\xi$")
        axs[1,0].set_ylabel(r"$r B \xi$")
        axs[1,0].set_xlabel(r"$\psi$")
        axs[1,1].set_xlabel(r"$z$")

        for a in np.ndarray.flatten(axs):
            a.axhline(0,color='C1',ls='--',lw=0.7)


        fig.suptitle(f"{self.filename}: t = {t}")
        #plt.contourf(zax, pax, xro.T, 20, cmap='plasma')
        #plt.colorbar()

    def plotXRO2D(self, t = -1):
        plt.figure()
        zax = self.zax[1:-1, 1:-1]
        pax = self.pax[1:-1, 1:-1]
        xro = self.xro2d[t]
        plt.contourf(zax, pax, xro.T, 20, cmap='plasma')
        plt.title(f"{self.filename} t = {t}")
        plt.colorbar()

    def plotGrowth2D(self, t=-1):
        fig, axs = plt.subplots(1,3,figsize=(14,4.5))
        zax = self.zax[1:-1, 1:-1]
        pax = self.pax[1:-1, 1:-1]
        ratio = self.xro2d[t] / self.xro2d[t-1]
        gamma = np.log(self.xro2d[t] / self.xro2d[t-1]) / self.dt

        def plot2d(f,ax,cmap='bwr', N=20, factor=1):
            mag = np.max(np.abs(f)) / factor
            levels = np.linspace(-mag,mag,N+1)
            c = ax.contourf(zax, pax, f.T, levels=levels, cmap=cmap, vmin=-mag, vmax=mag, extend='both')
            fig.colorbar(c, ax=ax, orientation='vertical')

        plot2d(self.xro2d[t], axs[0], factor=5)
        plot2d(ratio, axs[1], factor=20)
        plot2d(gamma, axs[2])

        axs[0].set_title(r"$\xi_t$")
        axs[1].set_title(r"$\xi_t/\xi_{t-1}$")
        axs[2].set_title(r"$\Delta t^{-1} \log(\xi_t/\xi_{t-1})$")
        axs[0].set_ylabel(r"$\psi$ index")
        axs[0].set_xlabel(r"$z$ index")

        fig.suptitle(f"{self.filename} : t = {t}")
        fig.tight_layout()
#        plt.colorbar()

    def plotEnergy(self):
        '''
        Plot energy conservation metrics
        '''

        fig,axs = plt.subplots(3,1, figsize=(10,8), sharex=True)
        axs[0].plot(self.tenergy)
        axs[1].plot(self.enkin)
        axs[2].plot(self.enpot)

        for a in axs:
            a.grid()
            a.set_yscale('symlog')

        axs[0].set_title("tenergy")
        axs[1].set_title("enkin")
        axs[2].set_title("enpot")
        axs[2].set_xlabel("time index")
        fig.suptitle(self.filename)

    def plotXit(self,iz=130,jp=1):
        '''
        Pick a single (z,psi) point and plot xi(t)
        '''
        fig,axs = plt.subplots(1,1)
        axs.plot(self.xro2d[:,iz,jp])
        axs.set_yscale('symlog')
        axs.set_title(f"{self.filename} (iz,jp) = {iz},{jp}")
        axs.set_xlabel("time step")
        axs.set_ylabel("xi")

    def plotBoundary(self):

        zax = self.z_axis[1:-1]
        pax = self.psi_axis[1:-1]
        fig,axs = plt.subplots(2,2, figsize=(10,8))
        f2,a = plt.subplots(2,2, figsize=(10,8))
        for t in [0,1,2,5,10,49,50]:
            x = self.xro2d[t]
            axs[0,0].plot(zax, self.xro2d[t,:,-1], label=t)
            axs[1,0].plot(zax, self.xro2d[t,:,0], label=t)
            axs[0,1].plot(pax, self.xro2d[t,-1,:], label=t)
            axs[1,1].plot(pax, self.xro2d[t,0,:], label=t)

            a[0,0].plot(zax, x[:,-1] - x[:,-2], label=t)
            a[1,0].plot(zax, x[:,1] - x[:,0], label=t)
            a[0,1].plot(pax, x[-1,:] - x[-1,:], label=t)
            a[1,1].plot(pax, x[1,:] - x[0,:], label=t)
        axs[0,0].set_title("Psi = -2")
        axs[1,0].set_title("Psi = 1")
        axs[0,1].set_title("Z = -2")
        axs[1,1].set_title("Z = 1")

        axs[0,0].legend()
        fig.suptitle(self.filename)
        fig.tight_layout()
