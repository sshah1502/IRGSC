from IRGSC._sgc import star_galaxy_classification
import numpy as np


def extinction_corrected_photometry(self):
    print("########################################################")
    print('Correcting the optical photometry of the probable stellar sources \
          for reddening and extinction.')
    print("")
    ps_phot = self.read_optical_data()
    ps1_objid, ps_ra, e_ps_ra, ps_dec, e_ps_dec, gpsf, e_gpsf, gkron, e_gkron, rpsf, e_rpsf, rkron, e_rkron, \
        ipsf, e_ipsf, ikron, e_ikron, zpsf, e_zpsf, zkron, e_zkron, ypsf, e_ypsf, ykron, e_ykron,\
            objinfoflag, qualityflag, ndetections, nstackdetections, ginfoflag, ginfoflag2, \
                ginfoflag3, rinfoflag, rinfoflag2, rinfoflag3,  iinfoflag, iinfoflag2, iinfoflag3,\
                    zinfoflag, zinfoflag2, zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3 = ps_phot
    
    e_gr = np.sqrt(e_gpsf**2 + e_rpsf**2)
    e_gi = np.sqrt(e_gpsf**2 + e_ipsf**2)
    e_gz = np.sqrt(e_gpsf**2 + e_zpsf**2)
    e_gy = np.sqrt(e_gpsf**2 + e_ypsf**2)
    e_ri = np.sqrt(e_rpsf**2 + e_ipsf**2)
    e_rz = np.sqrt(e_rpsf**2 + e_zpsf**2)
    e_ry = np.sqrt(e_rpsf**2 + e_ypsf**2)
    e_iy = np.sqrt(e_ipsf**2 + e_ypsf**2)
    e_iz = np.sqrt(e_ipsf**2 + e_zpsf**2)
    e_zy = np.sqrt(e_zpsf**2 + e_ypsf**2)

    """
    extinction in ps1 filters taken from Tonry et.al. 2012
    """
    ag = (self.ebv)*0.88*(3.613 - 0.0972*(gpsf - ipsf) + 0.01*(gpsf - ipsf)**2)
    ar = (self.ebv)*0.88*(2.585 - 0.0315*(gpsf - ipsf))
    ai = (self.ebv)*0.88*(1.908 - 0.0152*(gpsf - ipsf))
    az = (self.ebv)*0.88*(1.499 - 0.0023*(gpsf - ipsf))
    ay = (self.ebv)*0.88*(1.251 - 0.0027*(gpsf - ipsf))

    """
    error in extinction
    """
    e_ag = ((self.err_ebv)*(3.613 - 0.0972*(gpsf - ipsf) + 0.02*(gpsf - ipsf)**2))+\
            ((self.ebv)*((-0.0972*e_gi)+(0.02*(gpsf - ipsf)*e_gi)))
    e_ar = (self.err_ebv)*(2.585 - 0.0315*(gpsf - ipsf)) + (self.ebv)*(-0.0315*e_gi)
    e_ai = (self.err_ebv)*(1.908 - 0.0152*(gpsf - ipsf)) + (self.ebv)*(-0.0152*e_gi)
    e_az = (self.err_ebv)*(1.499 - 0.0023*(gpsf - ipsf)) + (self.ebv)*(-0.0023*e_gi)
    e_ay = (self.err_ebv)*(1.251 - 0.0027*(gpsf - ipsf)) + (self.ebv)*(-0.0027*e_gi)

    """
    prefix ec_ stands for extinction corrected magnitudes and e_ec_ stands 
    for error in those magnitudes
    """
    ec_gmag = gpsf - ag
    ec_rmag = ps1r - ar
    ec_imag = ipsf - ai
    ec_zmag = ps1z - az
    ec_ymag = ps1y - ay

    e_ec_gmag = np.sqrt((e_gpsf)**2 + (e_ag)**2)
    e_ec_rmag = np.sqrt((e_rpsf)**2 + (e_ar)**2)
    e_ec_imag = np.sqrt((e_ipsf)**2 + (e_ai)**2)
    e_ec_zmag = np.sqrt((e_zpsf)**2 + (e_az)**2)
    e_ec_ymag = np.sqrt((e_ypsf)**2 + (e_ay)**2)

    e_ec_gr = np.sqrt(e_ec_gmag**2 + e_ec_rmag**2)
    e_ec_gi = np.sqrt(e_ec_gmag**2 + e_ec_imag**2)
    e_ec_gz = np.sqrt(e_ec_gmag**2 + e_ec_zmag**2)
    e_ec_gy = np.sqrt(e_ec_gmag**2 + e_ec_ymag**2)
    e_ec_ri = np.sqrt(e_ec_rmag**2 + e_ec_imag**2)
    e_ec_rz = np.sqrt(e_ec_rmag**2 + e_ec_zmag**2)
    e_ec_ry = np.sqrt(e_ec_rmag**2 + e_ec_ymag**2)
    e_ec_iz = np.sqrt(e_ec_imag**2 + e_ec_zmag**2)
    e_ec_iy = np.sqrt(e_ec_imag**2 + e_ec_ymag**2)
    e_ec_zy = np.sqrt(e_ec_zmag**2 + e_ec_ymag**2)


    de_reddened_gr = ec_gmag - ec_rmag
    de_reddened_ri = ec_rmag - ec_imag
    de_reddened_gi = ec_gmag - ec_imag
    de_reddened_gy = ec_gmag - ec_ymag
    de_reddened_gz = ec_gmag - ec_zmag
    de_reddened_ry = ec_rmag - ec_ymag
    de_reddened_rz = ec_rmag - ec_zmag
    de_reddened_iy = ec_imag - ec_ymag
    de_reddened_iz = ec_imag - ec_zmag
    de_reddened_zy = ec_zmag - ec_ymag


    psf_phot = ps1_objid, ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, e_ec_gmag, gkron, e_gkron, \
        ec_rmag, e_ec_gmag, gkron, e_gkron, ec_imag, e_ec_gmag, gkron, e_gkron, ec_zmag,\
                e_ec_gmag, gkron, e_gkron, ec_ymag, e_ec_gmag, gkron, e_gkron, de_reddened_gr, \
                        de_reddened_gi, de_reddened_ri, de_reddened_gy, de_reddened_gz, \
                                de_reddened_ry, de_reddened_rz, de_reddened_iy, de_reddened_iz, \
                                        de_reddened_zy, e_ec_gr, e_ec_gi, e_ec_gz, e_ec_gy, e_ec_ri,\
                                                e_ec_rz, e_ec_ry, e_ec_iz, e_ec_iy, e_ec_zy, objinfoflag,\
                                                qualityflag, ndetections, nstackdetections, ginfoflag,\
                                                ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2, rinfoflag3,\
                                                iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2,\
                                                        zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3
    return psf_phot
