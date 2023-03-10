from .read_data import read_optical_data, read_nir_data
from ._extinction_correction import extinction_corrected_photometry
from ._sgc import star_galaxy_classification
from ._sam import read_sam_file, select_kurucz_models, select_phoenix_models
from ._fitting import find_nearest, calc_sf, reduced_chi2_all_stars, reduced_chi2_leq_2,\
    computed_reduced_chi2, compute_nir2, compute_flux, calc_sf2
from ._validate import validation_ukidss, plot_validation_plots
import .run_config

class irgsc():
    print('Starting the Program')
    print("")

    def __init__(self, ra, dec, ebv, err_ebv, aj, ah, ak, use_optimal_method=True, validate=False):
        self.use_optimal_method=use_optimal_method
        self.validate = validate
        self.ra = ra
        self.dec = dec
        self.ebv = ebv
        self.err_ebv = err_ebv
        self.aj = aj
        self.ah = ah
        self.ak = ak

#version
#author
#credit