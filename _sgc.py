from IRGSC.read_data import read_optical_data
import numpy as np
from matplotlib import pyplot as plt



def star_galaxy_classification(self):
    print("")
    print('Seperating the probable stellar sources from the optical dataset which is filtered\
     for nan values')
    print("")
    ps_phot = self.read_optical_data()
    print("")
    print('Using psf-kron criteria to seperate stars and galaxies')
    print("")
    
    ps1_objid, ps_ra, e_ps_ra, ps_dec, e_ps_dec, gpsf, e_gpsf, gkron, e_gkron, rpsf, e_rpsf, rkron, e_rkron, \
        ipsf, e_ipsf, ikron, e_ikron, zpsf, e_zpsf, zkron, e_zkron, ypsf, e_ypsf, ykron, e_ykron,\
            objinfoflag, qualityflag, ndetections, nstackdetections, ginfoflag, ginfoflag2, \
                ginfoflag3, rinfoflag, rinfoflag2, rinfoflag3,  iinfoflag, iinfoflag2, iinfoflag3,\
                    zinfoflag, zinfoflag2, zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3 = ps_phot
    sgc_index = np.where((gpsf - gkron < 0.05) & (rpsf - rkron < 0.05) & (ipsf - ikron < 0.05) & \
                         (zpsf - zkron < 0.05) & (ypsf - ykron < 0.05))[0]
    galaxy_index = np.where((gpsf - gkron > 0.05) & (rpsf - rkron > 0.05) & (ipsf - ikron > 0.05) & \
                            (zpsf - zkron > 0.05) & (ypsf - ykron > 0.05))[0]

    print('Number of probable stellar sources =' + str(len(sgc_index)) + 'and number of extended \
          sources =' + str(len(ipsf) - len(sgc_index)))
    print("")
    print('Now plotting the (g-r) vs (r-i) CCD which shows stars in a locus and galaxies \
          as random respectively')
    print("")

    plt.clf()
    plt.figure(figsize=(10,10))
    plt.scatter((gpsf[sgc_index] - rpsf[sgc_index]), (rpsf[sgc_index] - ipsf[sgc_index]), \
                s=5, color= 'm', alpha = 0.3, label='stellar sources')
    plt.scatter((gpsf[galaxy_index] - rpsf[galaxy_index]), (rpsf[galaxy_index] - ipsf[galaxy_index]), \
                s=5, color = 'k', alpha = 0.3, label = 'extended sources')
    plt.xlabel('$(g-r)$')
    plt.ylabel('$(r-i)$')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig('ccd_stars_and_galaxies_seperated.png')
    plt.clf()

    print('Plotting the (ipsf-ikron) vs (ikron) scatter plot which shows stars and \
          galaxies as magenta and black points respectively')
    print("")

    plt.clf()
    plt.figure(figsize=(10,10))
    plt.scatter(ipsf[sgc_index], (ipsf[sgc_index] - ikron[sgc_index]), s=5, color='m', alpha = 0.3\
                label='stellar sources')
    plt.scatter(ipsf[galaxy_index], (ipsf[galaxy_index] - ikron[galaxy_index]), s=5, \
                color='k', alpha = 0.3, label='extended sources')
    plt.xlabel('$i_{psf}$')
    plt.ylabel('$i_{psf}-i_{kron}$')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig('psf_vs_kron_stars_and_galaxies_seperated.png')
    plt.clf()

    ps_phot = ps1_objid[sgc_index], ps_ra[sgc_index], e_ps_ra[sgc_index], ps_dec[sgc_index],\
        e_ps_dec[sgc_index], gpsf[sgc_index], e_gpsf[sgc_index], gkron[sgc_index], e_gkron[sgc_index],\
            rpsf[sgc_index], e_rpsf[sgc_index], rkron[sgc_index], e_rkron[sgc_index], ipsf[sgc_index],\
                e_ipsf[sgc_index], ikron[sgc_index], e_ikron[sgc_index], zpsf[sgc_index], e_zpsf[sgc_index],\
                    zkron[sgc_index], e_zkron[sgc_index], ypsf[sgc_index], e_ypsf[sgc_index], ykron[sgc_index],\
                        e_ykron[sgc_index], objinfoflag[sgc_index], qualityflag[sgc_index], ndetections[sgc_index],\
                            nstackdetections[sgc_index], ginfoflag[sgc_index], ginfoflag2[sgc_index],\
                                ginfoflag3[sgc_index], rinfoflag[sgc_index], rinfoflag2[sgc_index], rinfoflag3[sgc_index],\
                                    iinfoflag[sgc_index], iinfoflag2[sgc_index], iinfoflag3[sgc_index],\
                                        zinfoflag[sgc_index], zinfoflag2[sgc_index], zinfoflag3[sgc_index],\
                                            yinfoflag[sgc_index], yinfoflag2[sgc_index], yinfoflag3[sgc_index]
    
    print("#####################################################")
    print('Created an input optical catalogue of stellar sources')
    print("######################################################")

    return psf_phot
