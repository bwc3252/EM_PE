# -*- coding: utf-8 -*-
'''
Me2017 model
----------------------
Models adapted from code found in `gwemlightcurves <https://github.com/mcoughlin/gwemlightcurves/tree/master/gwemlightcurves/KNModels/io>`_
'''

from __future__ import print_function
import numpy as np
import scipy.interpolate
from scipy.interpolate import interp1d

from .model import model_base

class me2017(model_base):
    '''
    Implementation of Metzger model

    Parameters
    ----------
    weight : float
        Weight of the model
    '''

    def __init__(self, weight=1):
        name = 'me2017'
        #param_names = ['mej', 'vej', 'dist']
        param_names = ['mej', 'vej']
        bands = ['u', 'g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
        model_base.__init__(self, name, param_names, bands, weight)
        self.mAB = None
        self.tdays = None

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters and run differential equation for lightcurve
        model. Currently uses questionable fixed time step method to integrate
        from tmin to tmax, then stores the lightcurve values for interpolation.

        Parameters
        ----------
        params : dict
            Dictionary mapping parameter names to their values
        t_bounds : list
            [upper bound, lower bound] pair for time values
        '''
        mej = params['mej']
        vej = params['vej']
        #dist = params['dist']
        dist = 40
        dt = 0.05
        self.tdays, self.mAB = self._calc_lc(0.5, t_bounds[1], dt, mej, vej, dist)

    def evaluate(self, t, band):
        '''
        Evaluate model at specific time values using the current parameters by
        interpolating the lightcurve calculated in the set_params() method.

        Parameters
        ----------
        tvec_days : np.ndarray
            Time values
        band : string
            Band to evaluate
        '''
        band_ind = dict(zip(self.bands, range(len(self.bands)))) # map bands to indices
        index = band_ind[band]
        lc = self.mAB[index]
        mask = np.isfinite(lc)
        f = interp1d(self.tdays[mask], lc[mask], fill_value='extrapolate')
        return f(t), 0

    def _calc_lc(self, tini, tmax, dt, mej, vej, dist):

        ### for now, hard-code some things
        beta = 3.0
        kappa_r = 1.0

        # ** define constants **
        c = 3.0e10
        mp = 1.67e-24
        Msun = 2.0e33
        kb = 1.38e-16
        sigSB = 5.67e-5
        h = 6.63e-27
        arad = 7.56e-15
        Mpc = 3.08e24

        # ** define parameters **

        # fiducial redshift/distance
        #z = 0.01
        #D = 39.5*Mpc
        z = 0.00
        D = 1e-5*Mpc

        # define desired observer band wavelengths (nm)
        # u (0), b (1), v (2), r (3), i (4), z (5), y(6), j (7), k (8), l (9)
        #lambdaobs = np.array([365., 445., 551., 658., 806., 900., 1020., 1220., 2190., 3450.])

        # u (0) g (1) r (2) i (3) z (4) y (5) J (6) H (7) K (8)
        lambdaobs = np.array([354.3, 477.56, 612.95, 748.46, 865.78, 960.31, 1235.0, 1662.0, 2159.0])

        nuobs = c/(1.0e-7*lambdaobs)
        nuobs = nuobs/(1.0 + z)

        # total ejecta mass
        M0 = mej*Msun
        # minimum initial velocity
        v0 = vej*c
        # velocity index (M ~ v**-beta)
        #beta = 3.
        # initial thermal energy of bulk
        E0 = (M0)*(v0**2.0)/2.0
        # normalization of opacity of r-process matter (~10 for lanthanides, ~1 for non-lanthanides)
        #kappa_r = 10.
        # IGNORE PARAMETERS BELOW THIS LINE
        # mass cut of free neutrons
        Mn = 1.0e-8*Msun
        # electron fraction & initial neutron mass fraction in outermost layers
        Ye = 0.1
        Xn0max = 1.0-2.0*Ye
        # engine (0 = off, 1 = on)
        engine_switch = 0
        # BH (0 = magnetar, 1 = BH)
        BH_switch = 0
        ej = 0.1
        # magnetar period (in seconds) and magnetic field (G)
        P = 0.7e-3
        B = 1.0e15
        # magnetar collapse time (in units of initial spin-down times)
        tcollapse = 10000000.

        # ** define time array in seconds **
        #tprec = 10000
        #tmin = np.log(0.1)
        #tmax = np.log(1.0e6)
        #t = np.arange(tprec)*(tmax-tmin)/(tprec-1.0) + tmin
        #t = np.exp(t)
        #tdays = t/(3600.*24.)

        tdays = np.arange(tini,tmax+dt,dt)
        t = tdays*(3600.*24.)
        tprec = len(t)

        # ** define mass/velocity array of outer ejecta, comprised of half of mass **
        mmin = np.log(1.0e-8)
        mmax = np.log(M0/Msun)
        mprec = 300
        m = np.arange(mprec)*(mmax-mmin)/(mprec-1.0) + mmin
        m = np.exp(m)

        #vm(where(m gt 0.5*M0/Msun)) = v0
        #vm(where(m le 0.5*M0/Msun)) = v0*(m(where(m le 0.5*M0/Msun))/(0.5*M0/Msun))^(-1./beta)
        vm = v0*(m/(M0/Msun))**(-1./beta)
        vm[vm > c] = c

        # define thermalization efficiency from Barnes+16
        # 1e-2 Msun, 0.2 c
        ca3 = 1.3
        cb3 = 0.2
        cd3 = 1.1
        # 1e-3, 0.3 c
        ca2 = 8.2
        cb2 = 1.2
        cd2 = 1.52
        # 1e-2, 0.1 c
        ca = 0.56
        cb = 0.17
        cd = 0.74
        eth = 0.36*(np.exp(-ca*tdays) + np.log(1.0+2*cb*(tdays**(cd)))/(2*cb*tdays**(cd)))
        eth2 = 0.36*(np.exp(-ca2*tdays) + np.log(1.0+2*cb2*(tdays**(cd2)))/(2*cb2*tdays**(cd2)))
        eth3 = 0.36*(np.exp(-ca3*tdays) + np.log(1.0+2*cb3*(tdays**(cd3)))/(2*cb3*tdays**(cd3)))

        # ** calculate magnetar power **
        Rns = 12.e5
        # moment of inertia
        Ins = 1.3e20
        Ins = Ins*1.0e25
        # magnetic moment
        mu = B*(Rns**(3.0))
        # angular rotation rate
        omega = 2.0*np.pi/P
        # rotational energy
        Erot = 0.5*Ins*omega**(2.0)
        # maximum spin-down luminosity
        Lsd0 = mu**(2.0)*(omega**(4.0))/c**(3.0)
        tsd0 = Erot/Lsd0
        Lsd = Lsd0/(1.0 + t/tsd0)**(2.0)
        Lsd[t > tcollapse*tsd0] = 0.0
        Lsd = Lsd/1.0e20
        Lsd = Lsd/1.0e20
        Lsd2 = Lsd

        if BH_switch:
            #*** calculate BH fall-back power
            Lsd = 2.0e11*(ej/0.1)*(t/0.1)**(-5./3.)
        if not engine_switch:
            Lsd[:] = 0.0

        # ** define diffusive mass depth (assumed beta = 3) **
        Mdiff = (4.0*np.pi*(M0)**(1./3.)*(v0*c*t**2.)/(3.0*kappa_r))**(3./4.)
        Mdiff[Mdiff > M0] = M0
        Mdiff = Mdiff/Msun

        # ** define radioactive heating rates **
        # neutron and r-process mass fractions
        Xn0 = Xn0max*2*np.arctan((Mn/(m*Msun))**(1.0))/np.pi
        Xr = 1.0-Xn0

        # define arrays in mass layer and time
        Xn = np.zeros((mprec,tprec))
        edotn = np.zeros((mprec,tprec))
        edotr = np.zeros((mprec,tprec))
        edot = np.zeros((mprec,tprec))
        kappa = np.zeros((mprec,tprec))
        kappan = np.zeros((mprec,tprec))
        kappar = np.zeros((mprec,tprec))

        # define specific heating rates and opacity of each mass layer
        t0 = 1.3
        sig = 0.11

        tarray = np.tile(t,(mprec,1))
        Xn0array = np.tile(Xn0,(tprec,1)).T
        Xrarray = np.tile(Xr,(tprec,1)).T
        etharray = np.tile(eth,(mprec,1))
        Xn = Xn0array*np.exp(-tarray/900.)
        edotn = 3.2e14*Xn
        edotr = 4.0e18*Xrarray*(0.5 - (1./np.pi)*np.arctan((tarray-t0)/sig))**(1.3)*etharray
        edotr = 2.1e10*etharray*((tarray/(3600.*24.))**(-1.3))
        edot = edotn + edotr
        kappan = 0.4*(1.0-Xn-Xrarray)
        kappar = kappa_r*Xrarray
        kappa = kappan + kappar

        # define total r-process heating of inner layer
        Lr = M0*4.0e18*(0.5 - (1./np.pi)*np.arctan((t-t0)/sig))**(1.3)*eth
        Lr = Lr/1.0e20
        Lr = Lr/1.0e20

        # *** define arrays by mass layer/time arrays ***
        ene = np.zeros((mprec,tprec))
        lum = np.zeros((mprec,tprec))
        lumpdv = np.zeros((mprec,tprec))
        lumedot = np.zeros((mprec,tprec))
        tdiff  = np.zeros((mprec,tprec))
        tau = np.zeros((mprec,tprec))
        # properties of photosphere
        Rphoto = np.zeros((tprec,))
        vphoto = np.zeros((tprec,))
        mphoto = np.zeros((tprec,))
        kappaphoto = np.zeros((tprec,))

        # *** define arrays for total ejecta (1 zone = deepest layer) ***
        # thermal energy
        E = np.zeros((tprec,))
        # kinetic energy
        Ek = np.zeros((tprec,))
        # velocity
        v = np.zeros((tprec,))
        R = np.zeros((tprec,))
        taues = np.zeros((tprec,))
        Lrad = np.zeros((tprec,))
        temp = np.zeros((tprec,))
        # setting initial conditions
        E[0] = E0/1.0e20
        E[0] = E[0]/1.0e20
        Ek[0] = E0/1.0e20
        Ek[0] = Ek[0]/1.0e20
        v[0] = v0
        R[0] = t[0]*v[0]
        dt = t[1:]-t[:-1]
        dm = m[1:]-m[:-1]

        marray = np.tile(m,(tprec,1)).T
        dmarray = np.tile(dm,(tprec,1)).T

        for j in xrange(tprec-1):
            # one zone calculation
            temp[j] = 1.0e10*(3.0*E[j]/(arad*4.0*np.pi*R[j]**(3.0)))**(0.25)
            if (temp[j] > 4000.):
                kappaoz = kappa_r
            if (temp[j] < 4000.):
                kappaoz = kappa_r*(temp[j]/4000.)**(5.5)
            kappaoz = kappa_r
            LPdV = E[j]*v[j]/R[j]
            tdiff0 = 3.0*kappaoz*M0/(4.0*np.pi*c*v[j]*t[j])
            tlc0 = R[j]/c
            tdiff0 = tdiff0+tlc0
            Lrad[j] = E[j]/tdiff0
            Ek[j+1] = Ek[j] + LPdV*(dt[j])
            v[j+1] = 1.0e20*(2.0*Ek[j]/(M0))**(0.5)
            E[j+1] = (Lr[j] + Lsd[j]-LPdV-Lrad[j])*(dt[j]) + E[j]
            R[j+1] = v[j+1]*(dt[j]) + R[j]
            taues[j+1] = (M0)*0.4/(4.0*R[j+1]**(2.0))

            #templayer = (3.0*ene[:-1,j]*dm*Msun/(arad*4.0*np.pi*(t[j]*vm[:-1])**(3.0)))**(0.25)
            #kappa_correction = np.ones(templayer.shape)
            #kappa_correction[templayer > 4000.] = 1.0
            #kappa_correction[templayer < 4000.] = 1.0*(templayer[templayer < 4000.]/4000.)**(5.5)
            #kappa_correction[:] = 1.0
            kappa_correction = 1.0

            tdiff[:-1,j] = 0.08*kappa[:-1,j]*m[:-1]*Msun*3*kappa_correction/(vm[:-1]*c*t[j]*beta)
            tau[:-1,j] = m[:-1]*Msun*kappa[:-1,j]/(4.0*np.pi*(t[j]*vm[:-1])**(2.0))
            lum[:-1,j] = ene[:-1,j]/(tdiff[:-1,j] + t[j]*(vm[:-1]/c))
            ene[:-1,j+1] = (edot[:-1,j] - (ene[:-1,j]/t[j]) - lum[:-1,j])*(dt[j]) + ene[:-1,j]
            lum[:-1,j] = lum[:-1,j]*(dm)*Msun

            tau[mprec-1,j] = tau[mprec-2,j]

            # photosphere
            pig1 = np.argmin(np.abs(tdiff[:,j]-t[j]))
            pig = np.argmin(np.abs(tau[:,j]-1.0))
            vphoto[j] = vm[pig]
            Rphoto[j] = vphoto[j]*t[j]
            mphoto[j] = m[pig]
            kappaphoto[j] = kappa[pig,j]


        Ltotm = np.sum(lum,axis=0)
        Ltotm = Ltotm/1.0e20
        Ltotm = Ltotm/1.0e20

        if engine_switch:
            Ltot = Lrad
            Tobs = 1.0e10*(Ltot/(4.0*np.pi*(R)**(2.0)*sigSB))**(0.25)
            if not BH_switch:
                tlife = (Lsd/1.0e5)**(0.5)*(v/(0.3*c))**(0.5)*(t/(3600.*24.))**(-0.5)
                Ltot = Ltot/(1.0+tlife)
        if not engine_switch:
            Ltot = Ltotm
            Tobs = 1.0e10*(Ltot/(4.0*np.pi*(Rphoto)**(2.0)*sigSB))**(0.25)

        nuobsarray = np.tile(nuobs,(tprec,1)).T
        expo = np.exp(h*nuobsarray/(kb*Tobs))-1.0
        F = (2.0*np.pi*(h*nuobsarray)*((nuobsarray/c)**(2.0))/expo)*(Rphoto/D)*(Rphoto/D)

        mAB = -2.5*np.log10(F) - 48.6

        mAB += 5*(np.log10(dist*1e6) - 1)

        mask = np.isfinite(mAB)

        return tdays, mAB
