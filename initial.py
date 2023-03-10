import os
import sys
import numpy as np

class IRGSC():
    from read_data import read_optical_data, read_nir_data
    from _extinction_correction import extinction_corrected_photometry
    from _sgc import star_galaxy_classification
    from _sam import read_sam_file, select_kurucz_models, select_phoenix_models
    from _fitting import find_nearest, calc_sf, reduced_chi2_all_stars, reduced_chi2_leq_2,\
        computed_reduced_chi2, compute_nir2, compute_flux, calc_sf2
    from _validate import validation_ukidss, plot_validation_plots
    
    def __init__(self, optical_data=None):
        print('Starting the program')
        self.optical_data = optical_data
        self.validating_data = None
        self.use_sam = None
        self.validate = False
        self.use_kurucz = False
        self.use_phoenix = False
        self.validate = False
        self.ra = mandatory_parameters[0]
        self.dec = mandatory_parameters[1]
        self.ebv = mandatory_parameters[2]
        self.e_ebv = mandatory_parameters[3]
        self.aj = mandatory_parameters[4]
        self.ah = mandatory_parameters[5]
        self.ak = mandatory_parameters[6]

        #self.use_interpolated_kurucz_models_greater_than_4000K = use_interpolated_kurucz_models_greater_than_4000K
        self.use_interpolated_kurucz_models_greater_than_4000K = False
        self.use_interpolated_kurucz_models_lesser_than_4000K = False
        self.use_interpolated_phoenix_models_greater_than_4000K = False
        self.use_interpolated_phoenix_models_lesser_than_4000K = False
        self.use_reduced_chi2_leq_2 = False
        self.use_reduced_chi2_all_stars = False
        self.use_compute_nir_by_keeping_sf_and_reddening_free = False
        self.starting_guess = [-21, 0.001] #starting guess for compute_nir2() subroutine
