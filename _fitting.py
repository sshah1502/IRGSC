import numpy as np
from IRGSC._sam import read_sam_file, select_sam
from IRGSC._extinction_correction import extinction_corrected_photometry
from matplotlib import pyplot as plt
from datetime import datetime
str_current_datetime = str(current_datetime)
field_name = str(self.field_name)

header = ['ps1_objid', 'ps1_ra', 'ps1_ra_error', 'ps1_dec', 'ps1_dec_error', 'ps1_gpsf', 'ps1_gpsf_error', \
'ps1_rpsf', 'ps1_rpsf_error', 'ps1_ipsf', 'ps1_ipsf_error', 'ps1_zpsf', 'ps1_zpsf_error', 'ps1_ypsf',\
    'ps1_ypsf_error', 'teff', 'logg', 'feh', 'sam_g','sam_r','sam_i','sam_z','sam_y','sam_j','sam_h',\
        'sam_k', 'scale_factor', 'scale_factor_error', 'chi2', 'computed_j', 'computed_j_error', 'computed_h',\
            'computed_h_error', 'computed_k', 'computed_k_error', 'gaia_source_id', 'gaia_ra', 'gaia_ra_error',\
                'gaia_dec', 'gaia_dec_error', 'gaia_parallax', 'gaia_parallax_error', 'gaia_pm', 'gaia_pm_ra',\
                    'gaia_pm_ra_error', 'gaia_pm_dec', 'gaia_pm_dec_error', 'gaia_ruwe', 'objinfoflag',\
                        'qualityflag', 'ndetections', 'nstackdetections', 'ginfoflag', 'ginfoflag2',\
                            'ginfoflag3', 'rinfoflag', 'rinfoflag2', 'rinfoflag3','iinfoflag', 'iinfoflag2',\
                                'iinfoflag3','zinfoflag', 'zinfoflag2', 'zinfoflag3', 'yinfoflag', 'yinfoflag2',\
                                    'yinfoflag3', 'SAM Flag']


def stdv(sfavg, v1, v2, v3, v4, v5): # function to calculate standard deviation
    mu = sfavg
    n=5
    sig = (((mu - v1)**2 + (mu - v2)**2 + (mu-v3)**2 + (mu - v4)**2 + (mu - v5)**2)/n)**0.5
    return sig

def find_nearest(self, array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def calc_sf(self, j, om, e_om, sm, index_min_ang_seperation):

    ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag = om
    e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag = e_om
    sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = sm

    sf_mean = (1/5.0)*((ec_gmag[j] - sam_g[index_min_ang_seperation]) + (ec_rmag[j] - sam_r[index_min_ang_seperation])\
               + (ec_imag[j] - sam_i[index_min_ang_seperation]) + (ec_zmag[j] - sam_z[index_min_ang_seperation])\
                 + (ec_ymag[j] - sam_y[index_min_ang_seperation]))

    e_sf_mean = (1/5)*np.sqrt(e_ec_gmag[j]**2 + e_ec_rmag[j]**2 + e_ec_imag[j]**2\
                                + e_ec_zmag[j]**2  + e_ec_ymag[j]**2 )

    """0.91 is the conversion constant from J_AB to J_Vega
        1.39 is the conversion constant from H_AB to H_Vega
        1.85 is the conversion constant from K_AB to K_Vega
    """

    cj = sf_mean + self.aj + sam_j[index_min_ang_seperation] - 0.91
    ch = sf_mean + self.ah + sam_h[index_min_ang_seperation] - 1.39
    ck = sf_mean + self.ak + sam_k[index_min_ang_seperation] - 1.85

    e_cj = np.sqrt(e_sf_mean**2)# + (self.e_aj)**2)
    e_ch = np.sqrt(e_sf_mean**2)# + (self.e_ah)**2)
    e_ck = np.sqrt(e_sf_mean**2)# + (self.e_ak)**2)

    return sf_mean, e_sf_mean, cj, e_cj, ch, e_ch, ck, e_ck


def computed_dquad(self, j, oc, mc, e_oc):
    """
        This routine calculates the quantity d_quad which is the cumulative
            difference in the observed and the model colors.
    """
    obs_gr, obs_gi, obs_gz, obs_gy, obs_ri, obs_ry, obs_rz, obs_iz, obs_iy, obs_zy = oc
    sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz, sam_ry, sam_iz, sam_iy, sam_zy = mc
    e_obs_gr, e_obs_ri, e_obs_gi, e_obs_gz, e_obs_gy, e_obs_rz, e_obs_ry, e_obs_iz, \
        e_obs_iy, e_obs_zy = e_oc
    dev_gr = (obs_gr[j] - sam_gr); dev_gi = (obs_gi[j] - sam_gi); dev_gz = (obs_gz[j] - sam_gz)
    dev_gy = (obs_gy[j] - sam_gy); dev_ri = (obs_ri[j] - sam_ri); dev_rz = (obs_rz[j] - sam_rz)
    dev_ry = (obs_ry[j] - sam_ry); dev_iz = (obs_iz[j] - sam_iz); dev_iy = (obs_iy[j] - sam_iy)
    dev_zy = (obs_zy[j] - sam_zy)
    dquad = (dev_gr**2 + dev_gi**2 + dev_gz**2 + dev_gy**2 + dev_ri**2\
                    + dev_rz**2 + dev_ry**2 + dev_iz**2 + dev_iy**2 + dev_zy**2)
    return dquad, np.min(dquad), dev_gr, dev_gi, dev_gz,\
                dev_gy, dev_ri, dev_rz, dev_ry, dev_iz, dev_iy, dev_zy


def model_fitting(self):

    dvf = []; model_params=[]; cat_ps_ra = []; cat_e_ps_ra = [];\
        cat_ps_dec = []; cat_e_ps_dec = []; cat_ec_gmag = []; cat_e_ec_gmag = [];\
        cat_ec_rmag = []; cat_e_ec_rmag = []; cat_ec_imag = []; cat_e_ec_imag = [];\
            cat_ec_zmag = []; cat_e_ec_zmag = []; cat_ec_ymag = []; cat_e_ec_ymag = [];\
            cat_sf = []; cat_e_sf = []; cat_min_dquad_element = []; cat_sam_g = []; cat_sam_r = [];\
                cat_sam_i = []; cat_sam_z = []; cat_sam_y = []; cat_teff = []; cat_logg = [];\
                cat_feh = []; cat_computed_j = []; cat_e_computed_j = []; cat_computed_h = [];\
                    cat_e_computed_h = []; cat_computed_k = []; cat_e_computed_k = []

    if self.use_optimal_method is True:
            print("")
            print('Computing the NIR magnitudes for all stars by computing chi2r after matching the colours')
            print("")
            ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, e_ec_gmag, gkron, e_gkron,\
                ec_rmag, e_ec_rmag, rkron, e_rkron, ec_imag, e_ec_imag, ikron, e_ikron, ec_zmag,\
                e_ec_zmag, zkron, e_zkron, ec_ymag, e_ec_ymag, ykron, e_ykron, de_reddened_gr,\
                    de_reddened_gi, de_reddened_ri, de_reddened_gy, de_reddened_gz, de_reddened_ry,\
                    de_reddened_rz, de_reddened_iy, de_reddened_iz, de_reddened_zy, e_ec_gr,\
                        e_ec_gi, e_ec_gz, e_ec_gy, e_ec_ri, e_ec_rz, e_ec_ry, e_ec_iz, e_ec_iy,\
                        e_ec_zy, objinfoflag, qualityflag, ndetections, nstackdetections,\
                            ginfoflag, ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2, rinfoflag3,\
                            iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2, zinfoflag3,\
                                yinfoflag, yinfoflag2, yinfoflag3 = self.extinction_corrected_photometry()

            gaia_source_id, gaia_ra, gaia_ra_error, gaia_dec, gaia_dec_error, gaia_parallax,\
                gaia_parallax_error, gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec,\
                    gaia_pm_dec_error, gaia_ruwe = self.obtain_gaia_data(self.ra,\
                                                                        self.dec, self.search_radius)

            model_params_k0 = self.select_sam(teff_range=[4000,10000],\
                                                         logg_range=None, feh_range=None)            
            model_params_c1 = self.select_sam(teff_range=[2800,5000],\
                                                          logg_range=[3.0,5.5], feh_range=[-5.0,-1.5])
            model_params_c2 = self.select_sam(teff_range=[2800,4000],\
                                                          logg_range=[0.0,3.0], feh_range=[-0.5,1.5])

            teff_c1, logg_c1, feh_c1, sam_g_c1, sam_r_c1, sam_i_c1, sam_z_c1, sam_y_c1, sam_j_c1,\
                sam_h_c1, sam_k_c1 = model_params_c1
            teff_c2, logg_c2, feh_c2, sam_g_c2, sam_r_c2, sam_i_c2, sam_z_c2, sam_y_c2, sam_j_c2,\
                sam_h_c2, sam_k_c2 = model_params_c1
            teff_k0, logg_k0, feh_k0, sam_g_k0, sam_r_k0, sam_i_k0, sam_z_k0, sam_y_k0, sam_j_k0,\
                sam_h_k0, sam_k_k0 = model_params_k0
            
            teff = np.append(teff_c1, teff_c2); teff = np.append(teff,teff_k0);\
                logg = np.append(logg_c1, logg_c2); logg = np.append(logg,logg_k0);\
                feh = np.append(feh_c1, feh_c2); feh = np.append(feh,feh_k0);\
                    sam_g = np.append(sam_g_c1, sam_g_c2); sam_g = np.append(sam_g,sam_g_k0);\
                    sam_r = np.append(sam_r_c1, sam_r_c2); sam_r = np.append(sam_r,sam_r_k0);\
                        sam_i = np.append(sam_i_c1, sam_i_c2); sam_i = np.append(sam_i,sam_i_k0);\
                        sam_z = np.append(sam_z_c1, sam_z_c2); sam_z = np.append(sam_z,sam_z_k0);\
                            sam_y = np.append(sam_y_c1, sam_y_c2); sam_y = np.append(sam_y,sam_y_k0);\
                            sam_j = np.append(sam_j_c1, sam_j_c2); sam_j = np.append(sam_j,sam_j_k0);\
                                sam_h = np.append(sam_h_c1, sam_h_c2); sam_h = np.append(sam_h,sam_h_k0);\
                                sam_k = np.append(sam_k_c1, sam_k_c2); sam_k = np.append(sam_k,sam_k_k0)

            sam_gr = sam_g - sam_r; sam_ri = sam_r - sam_i; sam_gi = sam_g - sam_i; sam_gz = sam_g - sam_z;\
                sam_gy = sam_g - sam_y; sam_ry = sam_r - sam_y; sam_rz = sam_r - sam_z; sam_iz = sam_i - sam_z;\
                    sam_iy = sam_i - sam_y; sam_zy = sam_z - sam_y

            observed_colours = de_reddened_gr, de_reddened_gi, de_reddened_gz, de_reddened_gy,\
                de_reddened_ri, de_reddened_ry, de_reddened_rz, de_reddened_iz, de_reddened_iy, \
                    de_reddened_zy
            
            model_colours = sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz, sam_ry, sam_iz, sam_iy, sam_zy
            
            e_observed_colours = e_ec_gr, e_ec_ri, e_ec_gi, e_ec_gz, e_ec_gy, e_ec_rz, e_ec_ry,\
                e_ec_iz, e_ec_iy, e_ec_zy

            observed_optical_magnitudes = ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag
            e_observed_optical_magnitudes = e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag
            sam_magnitudes = sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k

            with open('IRGSC_' + str(field_name) + '_' + str_current_datetime +'.csv', \
                      'w', encoding='UTF8') as file1:
                  writer=csv.writer(file1)
                  writer.writerow(header)
                  for j in range(len(de_reddened_gr)):
                    dquad_arr, min_dquad, _, _, _, _, _, _, _, _, _, _ = \
                        computed_dquad(j, oc = observed_colours, mc = model_colours, \
                                       e_oc = e_observed_colours)
                    min_dquad_element = find_nearest(dquad_arr, min_dquad)
                    index_best_fit_sam = np.where(min_dquad_element == (dquad_arr))[0]
                    sf_avg, sigma_sf, computed_j, computed_j_error, computed_h, computed_h_error, \
                        computed_k, computed_k_error = calc_sf(j=j, om = observed_optical_magnitudes, \
                                                            e_om = e_observed_optical_magnitudes,\
                                                            sm = sam_magnitudes, index_min_ang_seperation = index_best_fit_sam)
                    gaia_angular_seperation = 3600*np.sqrt(((ps_ra[j] - gaia_ra)\
                                                            * np.cos(np.radians(ps_dec[j])))**2\
                                                                + (ps_dec[j] - gaia_dec)**2)
                    index_min_ang_seperation = np.where(gaia_angular_seperation<=1.0)[0]
                    if len(index_min_ang_seperation) > 1.0:
                            gaia_ang_seperation_selected = gaia_angular_seperation[indg]
                            min_gaia_ang_seperation = gaia_angular_seperation\
                                [np.where(np.min(gaia_ang_seperation_selected)\
                                          == gaia_angular_seperation)[0]]

                            index_min_ang_seperation = np.where(min_gaia_ang_seperation \
                                                                == gaia_angular_seperation)[0]

                            data = ps1_objid[j], ps_ra[j], err_ps_ra[j], ps_dec[j], \
                                err_ps_dec[j], ec_gmag[j], e_ec_gmag[j], ec_rmag[j], \
                                e_ec_rmag[j], ec_imag[j], e_ec_imag[j], ec_zmag[j], \
                                    e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j], teff[index_min_ang_seperation][0], \
                                    logg[index_min_ang_seperation][0], feh[index_min_ang_seperation][0], sam_g[index_min_ang_seperation][0],\
                                        sam_r[index_min_ang_seperation][0], sam_i[index_min_ang_seperation][0], sam_z[index_min_ang_seperation][0],\
                                        sam_y[index_min_ang_seperation][0], sam_j[index_min_ang_seperation][0], sam_h[index_min_ang_seperation][0],\
                                            sam_k[index_min_ang_seperation][0], sf_avg[0], sigma_sf[j], min_dquad_element, computed_j[0],\
                                            computed_j_error[0], computed_h[0], computed_h_error[0], computed_k[0], computed_k_error[0],\
                                                gaia_source_id[index_min_ang_seperation][0], gaia_ra[index_min_ang_seperation][0],\
                                                gaia_ra_error[index_min_ang_seperation][0], gaia_dec[index_min_ang_seperation][0],\
                                                    gaia_dec_error[index_min_ang_seperation][0], gaia_parallax[index_min_ang_seperation][0],\
                                                    gaia_parallax_error[index_min_ang_seperation][0], gaia_pm[index_min_ang_seperation][0],\
                                                        gaia_pm_ra[index_min_ang_seperation][0], gaia_pm_ra_error[index_min_ang_seperation][0],\
                                                        gaia_pm_dec[index_min_ang_seperation][0], gaia_pm_dec_error[index_min_ang_seperation][0],\
                                                            gaia_ruwe[index_min_ang_seperation][0], objinfoflag[j], qualityflag[j], ndetections[j],\
                                                            nstackdetections[j], ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                                                                rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j], zinfoflag2[j],\
                                                                zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j]
                            writer.writerow(data)
                    elif len(index_min_ang_seperation) == 0.0:
                            data = ps1_objid[j], ps_ra[j], err_ps_ra[j], ps_dec[j], \
                                err_ps_dec[j], ec_gmag[j], e_ec_gmag[j], ec_rmag[j], \
                                e_ec_rmag[j], ec_imag[j], e_ec_imag[j], ec_zmag[j], \
                                    e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j], teff[index_min_ang_seperation][0], \
                                    logg[index_min_ang_seperation][0], feh[index_min_ang_seperation][0], sam_g[index_min_ang_seperation][0],\
                                        sam_r[index_min_ang_seperation][0], sam_i[index_min_ang_seperation][0], sam_z[index_min_ang_seperation][0],\
                                        sam_y[index_min_ang_seperation][0], sam_j[index_min_ang_seperation][0], sam_h[index_min_ang_seperation][0],\
                                            sam_k[index_min_ang_seperation][0], sf_avg[0], sigma_sf[j], min_dquad_element, computed_j[0],\
                                            computed_j_error[0], computed_h[0], computed_h_error[0], computed_k[0], computed_k_error[0],\
                                                -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, objinfoflag[j], qualityflag[j], ndetections[j],\
                                                nstackdetections[j], ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                                                    rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j], zinfoflag2[j],\
                                                    zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j]
                            writer.writerow(data)
                    elif len(index_min_ang_seperation) == 1.0:
                            data = ps1_objid[j], ps_ra[j], err_ps_ra[j], ps_dec[j], \
                                err_ps_dec[j], ec_gmag[j], e_ec_gmag[j], ec_rmag[j], \
                                e_ec_rmag[j], ec_imag[j], e_ec_imag[j], ec_zmag[j], \
                                    e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j], teff[index_min_ang_seperation][0], \
                                    logg[index_min_ang_seperation][0], feh[index_min_ang_seperation][0], sam_g[index_min_ang_seperation][0],\
                                        sam_r[index_min_ang_seperation][0], sam_i[index_min_ang_seperation][0], sam_z[index_min_ang_seperation][0],\
                                        sam_y[index_min_ang_seperation][0], sam_j[index_min_ang_seperation][0], sam_h[index_min_ang_seperation][0],\
                                            sam_k[index_min_ang_seperation][0], sf_avg[0], sigma_sf[j], min_dquad_element, computed_j[0],\
                                            computed_j_error[0], computed_h[0], computed_h_error[0], computed_k[0], computed_k_error[0],\
                                                gaia_source_id[index_min_ang_seperation][0], gaia_ra[index_min_ang_seperation][0],\
                                                gaia_ra_error[index_min_ang_seperation][0], gaia_dec[index_min_ang_seperation][0],\
                                                    gaia_dec_error[index_min_ang_seperation][0], gaia_parallax[index_min_ang_seperation][0],\
                                                    gaia_parallax_error[index_min_ang_seperation][0], gaia_pm[index_min_ang_seperation][0],\
                                                        gaia_pm_ra[index_min_ang_seperation][0], gaia_pm_ra_error[index_min_ang_seperation][0],\
                                                        gaia_pm_dec[index_min_ang_seperation][0], gaia_pm_dec_error[index_min_ang_seperation][0],\
                                                            gaia_ruwe[index_min_ang_seperation][0], objinfoflag[j], qualityflag[j], ndetections[j],\
                                                            nstackdetections[j], ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                                                                rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j], zinfoflag2[j],\
                                                                zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j]
                            writer.write(data)