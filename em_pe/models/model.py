from __future__ import print_function

class model_base:
    '''
    Base "template" class for lightcurve models
    '''
    def __init__(self, name, param_names, bands):
        self.name = name
        self.param_names = param_names
        self.bands = bands
        self.params = None

    def set_params(self, params):
        '''
        Method to set the parameters
        '''
        self.params = params

    ### Functions that should be implemented by child classes.
    ### NOTE: These could be seen as (and maybe should be) abstract methods

    def evaluate(self, t):
        '''
        Method to evaluate model at specific time values using the current parameters.

        Should return a dictionary of data band names mapped to values.
        '''
        pass
