def mandatory_parameters():
    self.ra = sys.argv[1]
    self.dec = sys.argv[2]
    self.ebv = sys.argv[3]
    self.err_ebv = sys.argv[4]
    self.aj = sys.argv[5]
    self.ah = sys.argv[6]
    self.ak = sys.argv[7]


def todo():
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
