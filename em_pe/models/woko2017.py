# -*- coding: utf-8 -*-
'''
Woko2017 Model
----------------------
Models adapted from code found in `gwemlightcurves <https://github.com/mcoughlin/gwemlightcurves/tree/master/gwemlightcurves/KNModels/io>`_
'''
from __future__ import print_function
import numpy as np
from scipy.interpolate import splrep, splev

from .model import model_base

class woko2017(model_base):
    '''
    Implementation of lightcurve model found `here <https://arxiv.org/abs/1705.07084>`_

    Parameters
    ----------
    weight : float
        Weight of the model
    '''

    def __init__(self, weight=1):
        name = 'woko2017'
        param_names = ['mej', 'vej', 'dist']
        bands = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
        model_base.__init__(self, name, param_names, bands, weight)
        modelfile = 'Data/DZ2_mags_2017-03-20.dat'
        #modelfile = 'Data/gamA2_mags_2017-03-20.dat'
        self.data_out = np.loadtxt(modelfile)
        self.band_ind = dict(zip(self.bands, range(len(self.bands)))) # map bands to indices
        self.f_list = []
        self.params = None

        ### Do the major interpolations
        ### define constants
        mej0 = 0.013+0.005
        vej0 = 0.132+0.08
        kappa0 = 1.0
        mejconst = [-1.13,-1.01,-0.94,-0.94,-0.93,-0.93,-0.95,-0.99]
        vejconst = [-1.28,-1.60,-1.52,-1.56,-1.61,-1.61,-1.55,-1.33]
        kappaconst = [2.65,2.27,2.02,1.87,1.76,1.56,1.33,1.13]
        ### temporarily hardcode these
        kappa_r = 1.0
        theta_r = 0.0
        data_out = self.data_out
        ndata, nslices = data_out.shape
        ints = np.arange(0,ndata,ndata/9)

        #mAB = np.zeros((len(tvec_days),8))

        for ii in xrange(len(ints)):
            idx = np.arange(ndata/9) + ii*(ndata/9)
            data_out_slice = data_out[idx,:]

            t = data_out_slice[:,1]
            data = data_out_slice[:,2:]

            nt, nbins = data.shape

            a_i = (360/(2*np.pi))*np.arccos(1 - np.arange(nbins)*2/float(nbins))
            b_i = (360/(2*np.pi))*np.arccos(1 - (np.arange(nbins)+1)*2/float(nbins))
            bins = (a_i + b_i)/2.0

            idx = np.argsort(np.abs(bins-theta_r*2*np.pi))
            idx1 = idx[0]
            idx2 = idx[1]
            weight1 = 1/np.abs(bins[idx1]-theta_r*2*np.pi)
            weight2 = 1/np.abs(bins[idx1]-theta_r*2*np.pi)
            if not np.isfinite(weight1):
                weight1, weight2 = 1.0, 0.0
            elif not np.isfinite(weight2):
                weight1, weight2 = 0.0, 1.0
            else:
                weight1, weight2 = weight1/(weight1+weight2), weight2/(weight1+weight2)

            if ii == 0:
                f1 = splrep(t, np.log10(data[:,idx1]))
                f2 = splrep(t, np.log10(data[:,idx2]))
            else:
                f1 = splrep(t,data[:,idx1])
                f2 = splrep(t,data[:,idx2])

            self.f_list.append((f1, f2, weight1, weight2))

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters for lightcurve
        model.

        Parameters
        ----------
        params : dict
            Dictionary mapping parameter names to their values
        t_bounds : list
            [upper bound, lower bound] pair for time values
        '''
        self.params = params


    def evaluate(self, tvec_days, band):
        '''
        Evaluate model at specific time values using the current parameters.

        Parameters
        ----------
        tvec_days : np.ndarray
            Time values
        band : string
            Band to evaluate
        '''
        ### define constants
        mej0 = 0.013+0.005
        vej0 = 0.132+0.08
        ### temporarily hardcode these
        kappa_r = 1.0
        theta_r = 0.0
        kappa0 = 1.0
        mejconst = [-1.13,-1.01,-0.94,-0.94,-0.93,-0.93,-0.95,-0.99]
        vejconst = [-1.28,-1.60,-1.52,-1.56,-1.61,-1.61,-1.55,-1.33]
        kappaconst = [2.65,2.27,2.02,1.87,1.76,1.56,1.33,1.13]

        ndata, nslices = self.data_out.shape

        mej = self.params['mej']
        vej = self.params['vej']
        #delta_t = self.params['delta_t']
        delta_t = 0

        #dist = 40.0

        # time shift
        tvec_days += delta_t

        ints = np.arange(0,ndata,ndata/9)
        mAB = np.zeros((len(tvec_days),8))

        for ii in xrange(len(ints)):
            f1, f2, weight1, weight2 = self.f_list[ii]
            fam1, fam2  = splev(tvec_days, f1), splev(tvec_days, f2)
            fam = weight1*fam1+weight2*fam2

            if ii == 0:
                lbol = 10**fam
            else:
                mAB[:,int(ii-1)] = np.squeeze(fam + mejconst[int(ii-1)]*np.log10(mej/mej0) + vejconst[int(ii-1)]*np.log10(vej/vej0))

        tmax = (kappa_r/10)**0.35 * (mej/10**-2)**0.318 * (vej/0.1)**-0.60
        Lmax = 2.8*10**40 * (kappa_r/10)**-0.60 * (mej/10**-2)**0.426 * (vej/0.1)**0.776

        tvec_days = tvec_days*tmax/tvec_days[np.argmax(lbol)]
        lbol = lbol*Lmax/np.max(lbol)
        ### u      g       r       i       z       y       J       H       K
        ### 354.3, 477.56, 612.95, 748.46, 865.78, 960.31, 1235.0, 1662.0, 2159.0
        wavelengths = [4775.6, 6129.5, 7484.6, 8657.8, 9603.1, 12350, 16620, 21590]
        wavelength_interp = 3543

        mAB_y = np.zeros(tvec_days.shape)
        for ii in xrange(len(tvec_days)):
            mAB_y[ii] = np.interp(wavelength_interp,wavelengths,mAB[ii,:])
        mAB_new = np.zeros((len(tvec_days),9))
        mAB_new[:,0] = np.squeeze(mAB_y)
        mAB_new[:,1:] = mAB
        #mAB_new += 5*(np.log10(dist*1e6) - 1)

        band_ind = dict(zip(self.bands, range(len(self.bands)))) # map bands to indices
        index = self.band_ind[band]
        return mAB_new.T[index], 0
