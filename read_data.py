import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

print("###########################")
print('Reading input datafiles')
print ("##########################")

def read_optical_data(self):

    print("########################")
    print('Reading optical data')
    print("########################")

    p1 = self.optical_data
    ps1_objid = p1[:,0]; ps_ra = p1[:,1]; ps_dec = p1[:,2]; ps_ra_err = p1[:,3]; ps_dec_err = p1[:,4];\
        gmag = p1[:,5]; gkron = p1[:,7]; e_gmag = p1[:,6]; e_gkron = p1[:,8]; rmag = p1[:,9];\
            rkron = p1[:,11]; e_rmag = p1[:,10]; e_rkron = p1[:,12]; imag = p1[:,13]; \
                ikron = p1[:,15]; e_imag = p1[:,14]; e_ikron = p1[:,16]; zmag = p1[:,17];\
                    zkron = p1[:,19]; e_zmag = p1[:,18]; e_zkron = p1[:,20]; ymag = p1[:,21];\
                        ykron = p1[:,23]; e_ymag = p1[:,22]; e_ykron = p1[:,24]; objinfoflag = p1[:,25];\
                            qualityflag = p1[:,26]; ndetections = p1[:,27]; nstackdetections = p1[:,28];\
                                ginfoflag = p1[:,29]; ginfoflag2 = p1[:,30]; ginfoflag3 = p1[:,31];\
                                    rinfoflag = p1[:,29]; rinfoflag2 = p1[:,30]; rinfoflag3 = p1[:,31];\
                                        iinfoflag = p1[:,29]; iinfoflag2 = p1[:,30]; iinfoflag3 = p1[:,31];\
                                            zinfoflag = p1[:,29]; zinfoflag2 = p1[:,30]; zinfoflag3 = p1[:,31];\
                                                yinfoflag = p1[:,29]; yinfoflag2 = p1[:,30]; yinfoflag3 = p1[:,31]

    print("")
    print('Now filtering the optical data for nan values')
    print("")

    ind = np.where((gpsf!= -999) & (ipsf!= -999) & (rpsf != -999) & (zpsf != -999) \
                   & (ypsf != -999) & (e_gpsf != -999) & (e_rpsf != -999) & (e_ipsf != -999) \
                    & (e_zpsf != -999) & (e_ypsf != -999) & (gkron!= -999) & (ikron!= -999) \
                        & (zkron != -999) & (ykron != -999) & (rkron != -999) & (e_gpsf < 0.2) \
                            & (e_rpsf < 0.2) & (e_ipsf < 0.2) & (e_zpsf < 0.2) & (e_ypsf < 0.2))[0]

    print('Number of Sources in the Optical Dataset=', len(ind))
    print("")

    print('Plotting the psf-kron vs psf scatter plot of the sources in the optical data')
    print("")

    plt.clf()
    plt.figure(figsize=(10,10))
    plt.scatter(ipsf[ind], ipsf[ind] - ikron[ind], color='green', s=5, alpha=0.5)
    plt.xlabel('$i_{psf}$')
    plt.ylabel('$i_{kron}$')
    plt.grid()
    plt.xlim(10,24)
    plt.ylim(-0.8,2.0)
    plt.savefig('ps1_all_sources_psf_vs_kron.png')
    plt.clf()

    print('Plotting the CCD of the sources in the optical data')
    print("")

    plt.clf()
    plt.scatter((gpsf[ind] - rpsf[ind]), (rpsf[ind] - ipsf[ind]), s=5,color='green', alpha = 0.3, label='Number of sources='+str(len(ind)))
    plt.xlabel('$(g-r)$')
    plt.ylabel('$(r-i)$')
    plt.grid()
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.2,2.2)
    plt.savefig('ps1_all_sources_ccd_whole.png')
    plt.clf()

    raw_optical_data = ps1_objid, ps_ra, ps_ra_err, ps_dec, ps_dec_err, gmag, e_gmag, gkron, e_gkron,\
                          rmag, e_rmag, rkron, e_rkron, imag, e_imag, ikron, e_ikron, zmag, e_zmag, zkron,\
                            e_zkron, ymag, e_ymag, ykron, e_ykron, objinfoflag, qualityflag, ndetections,\
                                nstackdetections, ginfoflag, ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2,\
                                    rinfoflag3,  iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2, \
                                        zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3
    return raw_optical_data

def read_nir_data(self):
    print('vd=', self.validate)
    if self.validate == True:
        if self.validating_data == None:
            print('Error: Validating data not provided')
        print('Now reading the NIR survey data for validation')
        #hdulist = fits.open('/mnt/c/Users/sshah/Documents/itcc/irgsc/script/tf6.fits',  memmap=True)
        validating_data = self.validating_data
        hdulist = validating_data#fits.open(str(validating_data), memap=True)

        p8 = hdulist[1].data

        petro_j = p8['JPETROMAG']
        e_petro_j = p8['jPetroMagErr']
        petro_h = p8['HPETROMAG']
        e_petro_h = p8['hPetroMagErr']
        petro_k = p8['KPETROMAG']
        e_petro_k = p8['kPetroMagErr']
        t1_ra = (p8['RA'])*(180.0/np.pi)
        t1_dec = p8['DEC']*(180.0/np.pi)

        nir_index = np.where((petro_j != -9.99999500e+08) & (e_petro_j != -9.99999500e+08) & (petro_h != -9.99999500e+08) & (e_petro_h != -9.99999500e+08) & (petro_k != -9.99999500e+08) & (e_petro_k != -9.99999500e+08) & (e_petro_j < 0.1) & (e_petro_h < 0.1) & (e_petro_k < 0.1))[0]
        print('Number of Stars in the NIR data = ', len(nir_index))

        filtered_petro_j = petro_j[nir_index]
        e_petro_h = e_petro_h[nir_index]
        filtered_petro_h = petro_h[nir_index]
        e_petro_k = e_petro_k[nir_index]
        filtered_petro_k = petro_k[nir_index]
        e_petro_j = e_petro_j[nir_index]
        filtered_t1_ra = t1_ra[nir_index]
        filtered_t1_dec = t1_dec[nir_index]

        plt.clf()
        plt.hist(filtered_petro_j)
        plt.xlabel('$J_{petro}$')
        plt.ylabel('Bins')
        plt.savefig('hist_ukidss_J.png')
        plt.clf()

        plt.clf()
        plt.hist(filtered_petro_h)
        plt.xlabel('$H_{petro}$')
        plt.ylabel('Bins')
        plt.savefig('hist_ukidss_H.png')
        plt.clf()

        plt.clf()
        plt.hist(filtered_petro_k)
        plt.xlabel('$K_{petro}$')
        plt.ylabel('Bins')
        plt.savefig('hist_ukidss_K.png')
        plt.clf()

        nir_data =  filtered_petro_j, filtered_petro_h, filtered_petro_k, e_petro_j, e_petro_h, e_petro_k, filtered_t1_ra, filtered_t1_dec
        return nir_data
