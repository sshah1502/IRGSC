import numpy as np

def read_sam_file(self, use_sam=None, teff_range=None, logg_range=None, feh_range=None):

        if self.use_sam == 'Kurucz':
            print("")
            print('Reading Interpolated Kurucz SAMs')
            print("")
            p2 = np.genfromtxt('interpolated_kurucz.txt')
        elif self.use_sam == 'Phoenix':
            print("")
            print('Reading Interpolated Phoenix SAMs')
            print("")
            p2 = np.genfromtxt('interpolated_phoenix.txt')
        
        teff = p2[:,0]; logg = p2[:,2]; feh = p2[:,1]; sam_g = p2[:,3]; sam_r = p2[:,4];\
                sam_i = p2[:,5]; sam_z = p2[:,6]; sam_y = p2[:,7]; sam_j = p2[:,8];\
                        sam_h = p2[:,9]; sam_k = p2[:,10]
        
        sam_params =  teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k
        return sam_params

def select_sam(self, teff_range=None, logg_range=None, feh_range=None):
        if self.use_optimal_method is True:
                if teff_range is None and logg_range is None and feh_range is None:
                        raise  TypeError("range must be provided in order to use the optimal method")
                elif teff_range is not None and logg_range is not None and feh_range is not None:
                        teff_lower_limit = teff_range[0]; teff_upper_limit = teff_range[-1];\
                        logg_lower_limit = logg_range[0]; logg_upper_limit = logg_range[-1];\
                        feh_lower_limit = feh_range[0]; feh_upper_limit = feh_range[-1]

                teffk, loggk, fehk, sam_gk, sam_rk, sam_ik, sam_zk, sam_yk, sam_jk, sam_hk, sam_kk\
                        = self.read_sam_file(use_sam='Kurucz')
                teffp, loggp, fehp, sam_gp, sam_rp, sam_ip, sam_zp, sam_yp, sam_jp, sam_hp, sam_kp\
                        = self.read_sam_file(use_sam='Phoenix')
                
                indk0 = np.where((teffk>teff_lower_limit) & (teffk<teff_upper_limit)\
                           & (fehk>feh_lower_limit) & (fehk<feh_upper_limit)\
                                & (loggk>logg_lower_limit) & (loggk<logg_upper_limit))

                teffk0, loggk0, fehk0, sam_gk0, sam_rk0, sam_ik0, sam_zk0, \
                        sam_yk0, sam_jk0, sam_hk0, sam_kk0 = teffk[indk0], loggk[indk0], \
                                fehk[indk0], sam_gk[indk0], sam_rk[indk0], sam_ik[indk0], \
                                        sam_zk[indk0], sam_yk[indk0], sam_jk[indk0], \
                                                sam_hk[indk0], sam_kk[indk0]
                
                indc0 = np.where((teffp>teff_lower_limit) & (teffp<teff_upper_limit)\
                           & (fehp>feh_lower_limit) & (fehp<feh_upper_limit)\
                                & (loggp>logg_lower_limit) & (loggp<logg_upper_limit))
                teffc0, loggc0, fehc0, sam_gc0, sam_rc0, sam_ic0, sam_zc0,\
                          sam_yc0, sam_jc0, sam_hc0, sam_kc0 = teffp[indc0], loggp[indc0], \
                                fehp[indc0], sam_gp[indc0], sam_rp[indc0], sam_ip[indc0], \
                                        sam_zp[indc0], sam_yp[indc0], sam_jp[indc0], \
                                                sam_hp[indc0], sam_kp[indc0]
                
                indc1 = np.where((teffp>teff_lower_limit) & (teffp<teff_upper_limit)\
                           & (fehp>feh_lower_limit) & (fehp<feh_upper_limit)\
                                & (loggp>logg_lower_limit) & (loggp<logg_upper_limit))
                teffc1, loggc1, fehc1, sam_gc1, sam_rc1, sam_ic1, sam_zc1,\
                          sam_yc1, sam_jc1, sam_hc1, sam_kc1 = teffp[indc1], loggp[indc1], \
                                fehp[indc1], sam_gp[indc1], sam_rp[indc1], sam_ip[indc1], \
                                        sam_zp[indc1], sam_yp[indc1], sam_jp[indc1], \
                                                sam_hp[indc1], sam_kp[indc1]\
                
                teff = np.append(teffk0, teffc0); teff = np.append(teff, teffc1)
                logg = np.append(loggk0, loggc0); logg = np.append(logg, loggc1)
                feh = np.append(fehk0, fehc0); feh = np.append(feh, fehc1)
                sam_g = np.append(sam_gk0, sam_gc0); sam_g = np.append(sam_g, sam_gc1)
                sam_r = np.append(sam_rk0, sam_rc0); sam_r = np.append(sam_r, sam_rc1)
                sam_i = np.append(sam_ik0, sam_ic0); sam_i = np.append(sam_i, sam_ic1)
                sam_z = np.append(sam_zk0, sam_zc0); sam_z = np.append(sam_z, sam_zc1)
                sam_y = np.append(sam_yk0, sam_yc0); sam_y = np.append(sam_y, sam_yc1)
                sam_j = np.append(sam_jk0, sam_jc0); sam_j = np.append(sam_j, sam_jc1)
                sam_h = np.append(sam_hk0, sam_hc0); sam_h = np.append(sam_h, sam_hc1)
                sam_k = np.append(sam_kk0, sam_kc0); sam_k = np.append(sam_k, sam_kc1)

                sam_params =  teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k
                return sam_params