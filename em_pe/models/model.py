from __future__ import print_function

class model_base:
    '''
    Base "template" class for lightcurve models

    name: string name of model

    param_names: [string] names of parameters

    bands: [string] names of data bands

    weight: float weight of model (defaults to 1, used as weight for multi-component likelihoods)
    '''
    def __init__(self, name, param_names, bands, weight=1):
        self.name = name
        self.param_names = param_names
        self.bands = bands
        self.weight = weight
        self.params = None
        self.t_bounds = None

    def set_params(self, params, t_bounds):
        '''
        Method to set the parameters. Would be overridden for more complex models,
        because this is where the differential equations would be solved.
        '''
        self.params = params
        self.t_bounds = t_bounds

    ### Functions that should be implemented by child classes.
    ### NOTE: These could be seen as (and maybe should be) abstract methods

    def evaluate(self, t, band):
        '''
        Method to evaluate model at specific time values using the current parameters.
        Must be implemented by all child classes.
        '''
        pass
