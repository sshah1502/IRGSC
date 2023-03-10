from _fitting import reduced_chi2_all_stars, reduced_chi2_leq_2
from read_data import read_nir_data
import numpy as np
from matplotlib import pyplot as plt


def validation_ukidss(self):
#Names of the columns:::::::::
#ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, ec_rmag, ec_imag,\
#ec_zmag, ec_ymag, e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag,\
#e_ec_ymag, teff, logg, feh, ukidss_j, e_ukidss_j, ukidss_h,\
#e_ukidss_h, ukidss_k, e_ukidss_k, computed_j, e_computed_j, computed_h,\
#e_computed_h, computed_k, e_computed_k

    print('Now validating the computed NIR magnitudes using the UKIDSS data in "tf6.fits" file')
    diff_jf = []; diff_hf = []; diff_kf = []; ob_j = []; ob_h = []; ob_k = []; ps_raf = []; ps_decf = []
    e_ps_raf = []; e_ps_decf = []; gmagf = []; rmagf = []; imagf = []; zmagf = []; ymagf = []
    tefff = []; loggf = []; fehf = []; e_gmagf = []; e_rmagf = []; e_imagf = []; e_zmagf = []
    e_ymagf = []; computed_jf = []; computed_hf = []; computed_kf = []; e_obs_jf = []; e_obs_hf = []
    e_obs_kf = []; e_computed_jf = []; e_computed_kf = []; e_computed_hf = []

    if self.validate is True:
        if self.validating_data is None:
            print('No validating data provided')
        else:
            ukidss_j, ukidss_h, ukidss_k, e_ukidss_j, e_ukidss_h, e_ukidss_k, ukidss_ra, ukidss_dec = self.read_nir_data()

        if self.use_reduced_chi2_leq_2 is True:
            ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, e_ec_gmag, ec_rmag, e_ec_rmag, ec_imag, e_ec_imag, ec_zmag, e_ec_zmag, ec_ymag, e_ec_ymag, sf, e_sf, minv, sam_g, sam_r, sam_i, sam_z, sam_y, teff, logg, feh, computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k = self.reduced_chi2_leq_2()
            print('ps_ra=', ps_ra)
        elif self.use_reduced_chi2_all_stars is True:
            ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, e_ec_gmag, ec_rmag, e_ec_rmag, ec_imag, e_ec_imag, ec_zmag , e_ec_zmag, ec_ymag, e_ec_ymag, sf, e_sf, minv, sam_g, sam_r, sam_i, sam_z, sam_y, teff, logg, feh, computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k = self.reduced_chi2_all_stars()
        elif self.use_compute_nir_by_keeping_sf_and_reddening_free is True:
            ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, e_ec_gmag, ec_rmag, e_ec_rmag, ec_imag, e_ec_imag, ec_zmag , e_ec_zmag, ec_ymag, e_ec_ymag, minv, sam_g, sam_r, sam_i, sam_z, sam_y, teff, logg, feh, computed_j, computed_h, computed_k = self.compute_nir2()
        else:
            print('Error: Catalogue conaining computed NIR magnitudes not found.')
        validate_params=[]
        with open('validated_catalogue.txt', 'w') as file4:
            for i1 in range(len(ps_ra)):
                #positionally matching the sources in the UKIDSS within 1" to the PS1 sources in the catalogue
                gamma = 3600*np.sqrt(((ps_ra[i1] - ukidss_ra)*np.cos(np.radians(ps_dec[i1])))**2 + (ps_dec[i1] - ukidss_dec)**2)
                index_position_match = np.where(gamma<=1.0)[0]
                if len(index_position_match) > 1:
                    matched_positions = gamma[index_position_match]
                    minimum_seperation = gamma[np.where(np.min(matched_positions) == gamma)[0]]
                    index_minimum_seperation = np.where(minimum_seperation == gamma)[0]

                    diff_j = ukidss_j[index_minimum_seperation] - computed_j[i1]
                    diff_h = ukidss_h[index_minimum_seperation] - computed_h[i1]
                    diff_k = ukidss_k[index_minimum_seperation] - computed_k[i1]

                    diff_jf = np.append(diff_jf, diff_j)
                    ob_j = np.append(ob_j, ukidss_j[index_minimum_seperation])
                    computed_jf = np.append(computed_jf, computed_j[i1])

                    diff_hf = np.append(diff_hf, diff_h)
                    ob_h = np.append(ob_h, ukidss_h[index_minimum_seperation])
                    computed_hf = np.append(computed_hf, computed_h[i1])

                    diff_kf = np.append(diff_kf, diff_k)
                    ob_k = np.append(ob_k, ukidss_k[index_minimum_seperation])
                    computed_kf = np.append(computed_kf, computed_k[i1])

                    ps_raf = np.append(ps_raf, ps_ra[i1])
                    ps_decf = np.append(ps_decf, ps_dec[i1])
                    e_ps_raf = np.append(e_ps_raf, e_ps_ra[i1])
                    e_ps_decf = np.append(e_ps_decf, e_ps_dec[i1])

                    gmagf = np.append(gmagf, ec_gmag[i1])
                    rmagf = np.append(rmagf, ec_rmag[i1])
                    imagf = np.append(imagf, ec_imag[i1])
                    zmagf = np.append(zmagf, ec_zmag[i1])
                    ymagf = np.append(ymagf, ec_ymag[i1])

                    e_gmagf = np.append(e_gmagf, e_ec_gmag[i1])
                    e_rmagf = np.append(e_rmagf, e_ec_rmag[i1])
                    e_imagf = np.append(e_imagf, e_ec_imag[i1])
                    e_zmagf = np.append(e_zmagf, e_ec_zmag[i1])
                    e_ymagf = np.append(e_ymagf, e_ec_ymag[i1])

                    tefff = np.append(tefff, teff[i1])
                    loggf = np.append(loggf, logg[i1])
                    fehf = np.append(fehf, feh[i1])

                    e_obs_jf = np.append(e_obs_jf, e_ukidss_j[index_minimum_seperation])
                    e_computed_jf = np.append(e_computed_jf, e_ukidss_j[index_minimum_seperation])

                    e_obs_hf = np.append(e_obs_hf, e_ukidss_h[index_minimum_seperation])
                    e_computed_hf = np.append(e_computed_hf, e_ukidss_h[index_minimum_seperation])

                    e_obs_kf = np.append(e_obs_kf, e_ukidss_k[index_minimum_seperation])
                    e_computed_kf = np.append(e_computed_kf, e_ukidss_k[index_minimum_seperation])
                    validate_params = ps_ra[i1], e_ps_ra[i1], ps_dec[i1], e_ps_dec[i1], ec_gmag[i1], ec_rmag[i1], ec_imag[i1], ec_zmag[i1], ec_ymag[i1], e_ec_gmag[i1], e_ec_rmag[i1], e_ec_imag[i1], e_ec_zmag[i1], e_ec_ymag[i1], teff[i1], logg[i1], feh[i1], ukidss_j[index_minimum_seperation], e_ukidss_j[index_minimum_seperation], ukidss_h[index_minimum_seperation], e_ukidss_h[index_minimum_seperation], ukidss_k[index_minimum_seperation], e_ukidss_k[index_minimum_seperation], computed_j[i1], e_computed_j[i1], computed_h[i1], e_computed_h[i1], computed_k[i1], e_computed_k[i1]
                    file4.write('%0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n'%(ps_ra[i1], e_ps_ra[i1], ps_dec[i1], e_ps_dec[i1], ec_gmag[i1], ec_rmag[i1], ec_imag[i1], ec_zmag[i1], ec_ymag[i1], e_ec_gmag[i1], e_ec_rmag[i1], e_ec_imag[i1], e_ec_zmag[i1], e_ec_ymag[i1], teff[i1], logg[i1], feh[i1], ukidss_j[index_minimum_seperation], e_ukidss_j[index_minimum_seperation], ukidss_h[index_minimum_seperation], e_ukidss_h[index_minimum_seperation], ukidss_k[index_minimum_seperation], e_ukidss_k[index_minimum_seperation], computed_j[i1], e_computed_j[i1], computed_h[i1], e_computed_h[i1], computed_k[i1], e_computed_k[i1]))

                elif len(index_position_match) == 1:
                    diff_j = ukidss_j[index_position_match] - computed_j[i1]
                    diff_h = ukidss_h[index_position_match] - computed_h[i1]
                    diff_k = ukidss_k[index_position_match] - computed_k[i1]

                    diff_jf = np.append(diff_jf, diff_j)
                    ob_j = np.append(ob_j, ukidss_j[index_position_match])
                    computed_jf = np.append(computed_jf, computed_j[i1])

                    diff_hf = np.append(diff_hf, diff_h)
                    ob_h = np.append(ob_h, ukidss_h[index_position_match])
                    computed_hf = np.append(computed_hf, computed_h[i1])

                    diff_kf = np.append(diff_kf, diff_k)
                    ob_k = np.append(ob_k, ukidss_k[index_position_match])
                    computed_kf = np.append(computed_kf, computed_k[i1])

                    ps_raf = np.append(ps_raf, ps_ra[i1])
                    ps_decf = np.append(ps_decf, ps_dec[i1])
                    e_ps_raf = np.append(e_ps_raf, e_ps_ra[i1])
                    e_ps_decf = np.append(e_ps_decf, e_ps_dec[i1])

                    gmagf = np.append(gmagf, ec_gmag[i1])
                    rmagf = np.append(rmagf, ec_rmag[i1])
                    imagf = np.append(imagf, ec_imag[i1])
                    zmagf = np.append(zmagf, ec_zmag[i1])
                    ymagf = np.append(ymagf, ec_ymag[i1])

                    e_gmagf = np.append(e_gmagf, e_ec_gmag[i1])
                    e_rmagf = np.append(e_rmagf, e_ec_rmag[i1])
                    e_imagf = np.append(e_imagf, e_ec_imag[i1])
                    e_zmagf = np.append(e_zmagf, e_ec_zmag[i1])
                    e_ymagf = np.append(e_ymagf, e_ec_ymag[i1])

                    tefff = np.append(tefff, teff[i1])
                    loggf = np.append(loggf, logg[i1])
                    fehf = np.append(fehf, feh[i1])

                    e_obs_jf = np.append(e_obs_jf, e_ukidss_j[index_position_match])
                    e_computed_jf = np.append(e_computed_jf, e_ukidss_j[index_position_match])

                    e_obs_hf = np.append(e_obs_hf, e_ukidss_h[index_position_match])
                    e_computed_hf = np.append(e_computed_hf, e_ukidss_h[index_position_match])

                    e_obs_kf = np.append(e_obs_kf, e_ukidss_k[index_position_match])
                    e_computed_kf = np.append(e_computed_kf, e_ukidss_k[index_position_match])
                    validate_params = ps_ra[i1], e_ps_ra[i1], ps_dec[i1], e_ps_dec[i1], ec_gmag[i1], ec_rmag[i1], ec_imag[i1], ec_zmag[i1], ec_ymag[i1], e_ec_gmag[i1], e_ec_rmag[i1], e_ec_imag[i1], e_ec_zmag[i1], e_ec_ymag[i1], teff[i1], logg[i1], feh[i1], ukidss_j[index_position_match], e_ukidss_j[index_position_match], ukidss_h[index_position_match], e_ukidss_h[index_position_match], ukidss_k[index_position_match], e_ukidss_k[index_position_match], computed_j[i1], e_computed_j[i1], computed_h[i1], e_computed_h[i1], computed_k[i1], e_computed_k[i1]
                    file4.write('%0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n'%(ps_ra[i1], e_ps_ra[i1], ps_dec[i1], e_ps_dec[i1], ec_gmag[i1], ec_rmag[i1], ec_imag[i1], ec_zmag[i1], ec_ymag[i1], e_ec_gmag[i1], e_ec_rmag[i1], e_ec_imag[i1], e_ec_zmag[i1], e_ec_ymag[i1], teff[i1], logg[i1], feh[i1], ukidss_j[index_position_match], e_ukidss_j[index_position_match], ukidss_h[index_position_match], e_ukidss_h[index_position_match], ukidss_k[index_position_match], e_ukidss_k[index_position_match], computed_j[i1], e_computed_j[i1], computed_h[i1], e_computed_h[i1], computed_k[i1], e_computed_k[i1]))
        return validate_params

import os.path

def plot_validation_plots(self):
    if self.validate is True:
        file_exists = os.path.exists('validated_catalogue.txt')
        print('validated catalogue file exists ?', file_exists)
        if file_exists is False:
            print('validated catalogue does not exist')
            print('Now running validation with observed NIR data')
            self.validation_ukidss()
        elif file_exists is True:
            print('Reading the already existing validated_catalogue.txt file')
            p=np.genfromtxt('validated_catalogue.txt')
            ps_ra = p[:,0]; e_ps_ra = p[:,1]; ps_dec = p[:,2]; e_ps_dec = p[:,3];\
            ec_gmag = p[:,4]; ec_rmag = p[:,5]; ec_imag = p[:,6]; ec_zmag = p[:,7];\
            ec_ymag = p[:,8]; e_ec_gmag = p[:,9]; e_ec_rmag = p[:,10]; e_ec_imag = p[:,11];\
            e_ec_zmag = p[:,12]; e_ec_ymag = p[:,13]; teff = p[:,14]; logg = p[:,15]; feh = p[:,16];\
            ob_j = p[:,17]; e_obs_jf = p[:,18]; ob_h = p[:,19]; e_obs_hf = p[:,20]; ob_k = p[:,21];\
            e_obs_kf = p[:,22]; computed_jf = p[:,23]; e_computed_jf = p[:,24];\
            computed_hf = p[:,25]; e_computed_hf = p[:,26]; computed_kf = p[:,27];\
            e_computed_kf = p[:,28]

            diff_jf = ob_j - computed_jf
            diff_hf = ob_h - computed_hf
            diff_kf = ob_k - computed_kf

            plt.clf()
            plt.scatter(ob_j, (e_computed_jf**2)**0.5, s=5, alpha = 0.5)
            plt.grid()
            plt.ylabel('Error in $J_{computed}$')
            plt.xlabel('$J_{UKIDSS}$')
            plt.savefig('obj_vs_err_computed_j.png')
            plt.clf()

            plt.clf()
            plt.scatter(ob_j, (e_obs_jf**2)**0.5, s=5, alpha = 0.5)
            plt.grid()
            plt.ylabel('Error in $J_{UKIDSS}$')
            plt.xlabel('$J_{UKIDSS}')
            plt.savefig('obj_vs_err_obj.png')
            plt.clf()

            plt.clf()
            plt.scatter(ob_h, (e_computed_hf**2)**0.5, s=5, alpha = 0.5)
            plt.grid()
            plt.ylabel('Error in $H_{computed}$')
            plt.xlabel('$H_{UKIDSS}$')
            plt.savefig('obh_vs_err_computed_h.png')
            plt.clf()

            plt.clf()
            plt.scatter(ob_h, (e_obs_hf**2)**0.5, s=5, alpha = 0.5)
            plt.grid()
            plt.ylabel('Error in $H_{UKIDSS}$')
            plt.xlabel('$H_{UKIDSS}')
            plt.savefig('obh_vs_err_obh.png')
            plt.clf()

            plt.clf()
            plt.scatter(ob_k, (e_computed_kf**2)**0.5, s=5, alpha = 0.5)
            plt.grid()
            plt.ylabel('Error in $K_{computed}$')
            plt.xlabel('$K_{UKIDSS}$')
            plt.savefig('obk_vs_err_computed_k.png')
            plt.clf()

            plt.clf()
            plt.scatter(ob_k, (e_obs_kf**2)**0.5, s=5, alpha = 0.5)
            plt.grid()
            plt.ylabel('Error in $K_{UKIDSS}$')
            plt.xlabel('$K_{UKIDSS}')
            plt.savefig('obk_vs_err_obk.png')
            plt.clf()


            bins2 = np.arange(diff_jf.min(), diff_jf.max()+.1, 0.2)
            plt.clf()
            fig = plt.figure(figsize=(8,8))
            from matplotlib.gridspec import GridSpec
            gs = GridSpec(4,4)
            ax_joint = fig.add_subplot(gs[1:4,0:3])
            ax_marg_x = fig.add_subplot(gs[0,0:3])
            ax_marg_y = fig.add_subplot(gs[1:4,3])
            ax_joint.scatter(ob_j, diff_jf, alpha = 0.3, s = 5, color = 'g', label = 'No. of stars =' + str(len(ob_j)))
            ax_joint.grid()
            ax_joint.set_ylim(-2,2)
            ax_joint.legend(fontsize=14, loc = 'best')
            nx, bx, px = ax_marg_x.hist(ob_j, color = 'm', edgecolor = 'g', alpha = 0.5, label = 'Observed J')
            ny, by, px = ax_marg_y.hist(diff_jf, bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Difference')
            biny_max = np.where(ny == ny.max())
            print('maxbin', "{:.2f}".format(by[biny_max][0]))
            plt.text(200, -1.5, str('mode at '"{:.2f}".format(by[biny_max][0])), rotation = 270)
            ax_marg_y.set_ylim(-2,2)
            ax_marg_x.grid()
            ax_marg_x.legend(loc='best')
            ax_marg_y.grid()
            ax_marg_y.legend(loc='best')
            # Turn off tick labels on marginals
            plt.setp(ax_marg_x.get_xticklabels(), visible=False)
            plt.setp(ax_marg_y.get_yticklabels(), visible=False)
            # Set labels on joint
            ax_joint.set_xlabel('Observed J magnitude')
            ax_joint.set_ylabel('$J_{UKIDSS}$ - $J_{Computed}$')
            # Set labels on marginals
            ax_marg_y.set_xlabel('N')
            ax_marg_x.set_ylabel('N')
            ax_marg_x.set_title('Validation plot $J_{Computed}$')
            plt.savefig('validation_plot_j.png')
            plt.clf()

            bins2 = np.arange(diff_hf.min(), diff_hf.max()+.1, 0.2)
            plt.clf()
            fig = plt.figure(figsize=(8,8))
            from matplotlib.gridspec import GridSpec
            gs = GridSpec(4,4)
            ax_joint = fig.add_subplot(gs[1:4,0:3])
            ax_marg_x = fig.add_subplot(gs[0,0:3])
            ax_marg_y = fig.add_subplot(gs[1:4,3])
            ax_joint.scatter(ob_h, diff_hf, alpha = 0.3, s = 5, color = 'g', label = 'No. of stars =' + str(len(ob_j)))
            ax_joint.grid()
            ax_joint.set_ylim(-2,2)
            ax_joint.legend(fontsize=14, loc = 'best')
            nx, bx, px = ax_marg_x.hist(ob_h, color = 'm', edgecolor = 'g', alpha = 0.5, label = 'Observed J')
            ny, by, px = ax_marg_y.hist(diff_hf, bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Difference')
            biny_max = np.where(ny == ny.max())
            print('maxbin', "{:.2f}".format(by[biny_max][0]))
            plt.text(200, -1.5, str('mode at '"{:.2f}".format(by[biny_max][0])), rotation = 270)
            ax_marg_y.set_ylim(-2,2)
            ax_marg_x.grid()
            ax_marg_x.legend(loc='best')
            ax_marg_y.grid()
            ax_marg_y.legend(loc='best')
            # Turn off tick labels on marginals
            plt.setp(ax_marg_x.get_xticklabels(), visible=False)
            plt.setp(ax_marg_y.get_yticklabels(), visible=False)
            # Set labels on joint
            ax_joint.set_xlabel('Observed H magnitude')
            ax_joint.set_ylabel('$H_{UKIDSS}$ - $H_{Computed}$')
            # Set labels on marginals
            ax_marg_y.set_xlabel('N')
            ax_marg_x.set_ylabel('N')
            plt.title('Validation plot $H_{Computed}$')
            plt.savefig('validation_plot_h.png')
            plt.clf()

            bins2 = np.arange(diff_kf.min(), diff_kf.max()+.1, 0.2)
            plt.clf()
            fig = plt.figure(figsize=(8,8))
            from matplotlib.gridspec import GridSpec
            gs = GridSpec(4,4)
            ax_joint = fig.add_subplot(gs[1:4,0:3])
            ax_marg_x = fig.add_subplot(gs[0,0:3])
            ax_marg_y = fig.add_subplot(gs[1:4,3])
            ax_joint.scatter(ob_k, diff_kf, alpha = 0.3, s = 5, color = 'g', label = 'No. of stars =' + str(len(ob_j)))
            ax_joint.grid()
            ax_joint.set_ylim(-2,2)
            ax_joint.legend(fontsize=14, loc = 'best')
            nx, bx, px = ax_marg_x.hist(ob_k, color = 'm', edgecolor = 'g', alpha = 0.5, label = 'Observed J')
            ny, by, px = ax_marg_y.hist(diff_kf, bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Difference')
            biny_max = np.where(ny == ny.max())
            print('maxbin', "{:.2f}".format(by[biny_max][0]))
            plt.text(200, -1.5, str('mode at '"{:.2f}".format(by[biny_max][0])), rotation = 270)
            ax_marg_y.set_ylim(-2,2)
            ax_marg_x.grid()
            ax_marg_x.legend(loc='best')
            ax_marg_y.grid()
            ax_marg_y.legend(loc='best')
            # Turn off tick labels on marginals
            plt.setp(ax_marg_x.get_xticklabels(), visible=False)
            plt.setp(ax_marg_y.get_yticklabels(), visible=False)
            # Set labels on joint
            ax_joint.set_xlabel('Observed K magnitude')
            ax_joint.set_ylabel('$K_{UKIDSS}$ - $K_{Computed}$')
            # Set labels on marginals
            ax_marg_y.set_xlabel('N')
            ax_marg_x.set_ylabel('N')
            plt.title('Validation plot $K_{Computed}$')
            plt.savefig('validation_plot_k.png')
            plt.clf()
