"""This module handles calculations for a single gravitational wave detector.

"""

from copy import copy

import dill
import numpy as np

import gwbench.antenna_pattern_np as ant_pat_np
import gwbench.basic_functions as bfs
import gwbench.detector_responses as dr
import gwbench.err_deriv_handling as edh
import gwbench.fisher_analysis_tools as fat
import gwbench.psd as psd
import gwbench.snr as snr_mod

class Detector:

    ###
    #-----Init methods-----
    def __init__(self,det_key):
        #-----detector specification-----
        # full detector specification, specifying technology and location, e.g. CE2-40-CBO_C
        self.det_key = det_key
        # detector technology and locations
        self.tec = det_key.split('_')[0]
        self.loc = det_key.split('_')[1]

        #-----waveform and injection based quantities-----
        # frequency array
        self.f = None

        #-----technology and location based quantities-----
        # detector PSD
        self.psd = None
        # antenna pattern
        self.Fp = None
        self.Fc = None
        # location phase factor
        self.Flp = None

        #-----detector reponses-----
        # detector repsonse
        self.hf = None
        # derivative dictionary for detector responses
        self.del_hf = None
        # sympy expression of derivative dictionary for detector responses
        self.del_hf_expr = None

        #-----SNR-----
        # SNR, SNR^2 and d(SNR^2) calculated from self.hf
        self.snr = None
        self.snr_sq = None
        self.d_snr_sq = None

        #-----errors-----
        # Fisher matrix
        self.fisher = None
        # condition number of Fisher matrix
        self.cond_num = None
        # well conditioned Fisher matrix (True for yes, False for no)
        self.wc_fisher = None
        # covariance matrix
        self.cov = None
        # inversion error between the two matrices
        self.inv_err = None
        # dictionary of errors for given derivative variables
        self.errs = None


    ###
    #-----Setter methods-----
    def set_f(self, f):
        self.f = copy(f)


    ###
    #-----PSDs and antenna patterns-----
    def setup_ant_pat_lpf_psds(self, inj_params, use_rot, F_lo=-np.inf, F_hi=np.inf, psd_file_dict=None):
        self.setup_psds(F_lo, F_hi, psd_file_dict)
        self.setup_ant_pat_lpf(inj_params, use_rot)

    def setup_psds(self, F_lo=-np.inf, F_hi=np.inf, psd_file_dict=None):
        if psd_file_dict is None or self.det_key not in list(psd_file_dict.keys()):
            psd_file = None
            is_asd   = None
        else:
            psd_file = psd_file_dict[self.det_key]['psd_file']
            is_asd   = psd_file_dict[self.det_key]['is_asd']
        self.psd, self.f = psd.psd(self.tec,self.f,F_lo,F_hi,psd_file,is_asd)

    def setup_ant_pat_lpf(self, inj_params, use_rot):
        if use_rot:
            self.Fp, self.Fc, self.Flp = ant_pat_np.antenna_pattern_and_loc_phase_fac(
                self.f, inj_params['Mc'], inj_params['tc'],
                inj_params['ra'], inj_params['dec'], inj_params['psi'],
                inj_params['gmst0'], self.loc, use_rot)
        else:
            self.Fp, self.Fc, self.Flp = ant_pat_np.antenna_pattern_and_loc_phase_fac(
                self.f, None, None,
                inj_params['ra'], inj_params['dec'], inj_params['psi'],
                inj_params['gmst0'], self.loc, use_rot)


    ###
    #-----Detector responses-----
    def calc_det_responses(self, wf, inj_params):
        hfp, hfc = wf.eval_np_func(self.f,bfs.get_sub_dict(inj_params,wf.wf_symbs_string))
        self.hf = self.Flp * (hfp * self.Fp + hfc * self.Fc)

    def calc_det_responses_derivs_num(self, inj_params, deriv_variables, wf, wf_deriv_symbs_string, conv_cos, conv_log, deriv_symbs_string, use_rot, step, method, order, n):
        print(' ',self.det_key)
        self.calc_det_responses(wf,inj_params)
        self.del_hf = dr.calc_det_responses_num(self.loc,wf,wf_deriv_symbs_string,self.f,inj_params,use_rot,'hf',step,method,order,n)
        self.del_hf, c_quants = get_conv_del_eval_dic(self.del_hf, inj_params, conv_cos, conv_log, deriv_symbs_string)
        inj_params, deriv_variables = get_conv_inj_params_deriv_variables(c_quants, inj_params, deriv_variables)

    def load_det_responses_derivs_sym(self, wf_model_name, deriv_symbs_string, return_bin=0):
        self.del_hf_expr = dr.load_det_response_sym(self.loc, wf_model_name, deriv_symbs_string, return_bin)

    def calc_det_responses_derivs_sym(self, wf, inj_params, deriv_variables, conv_cos, conv_log, deriv_symbs_string):
        print(' ',self.det_key)
        self.calc_det_responses(wf,inj_params)
        self.del_hf = {}
        for deriv in self.del_hf_expr:
            if deriv in ('variables','deriv_variables'): continue
            self.del_hf[deriv] = self.del_hf_expr[deriv](self.f, **bfs.get_sub_dict(inj_params, self.del_hf_expr['variables']))

        self.del_hf, c_quants = get_conv_del_eval_dic(self.del_hf, inj_params, conv_cos, conv_log, deriv_symbs_string)
        inj_params, deriv_variables = get_conv_inj_params_deriv_variables(c_quants, inj_params, deriv_variables)


    ###
    #-----SNR calculations-----
    def calc_snrs_det_responses(self, only_net, df):
        print(' ',self.det_key)
        snr,snr_sq = snr_mod.snr_snr_sq_freq_array(self.hf, self.psd, self.f, df)
        if not only_net:
            self.snr = snr
            self.snr_sq = snr_sq
        return snr_sq

    def calc_snr_sq_integrand_det_responses(self):
        print(' ',self.det_key)
        self.d_snr_sq = snr_mod.snr_square_integrand(self.hf, self.psd)


    ###
    #-----Error calculation and Fisher analysis-----
    def calc_error_mats(self, only_net, df, cond_sup):
        print(' ',self.det_key)
        del_vs_f_dic = bfs.get_sub_dict(self.del_hf,('hf',),0)
        if not only_net:
            self.fisher, self.cov, self.wc_fisher, self.cond_num = fat.fisher_cov_matrix(list(del_vs_f_dic.values()),self.psd,self.f,0,df,cond_sup)
            return self.fisher
        else:
            fisher,_,_,_ = fat.fisher_cov_matrix(list(del_vs_f_dic.values()),self.psd,self.f,1,df,cond_sup)
            return fisher

    def calc_inv_err(self,by_element=0):
        self.inv_err = fat.inv_err_from_fisher_cov(self.fisher,self.cov,by_element)

    def calc_errs(self, deriv_variables):
        self.errs = fat.get_errs_from_cov(self.cov,deriv_variables)

    def calc_sky_area_90(self, ra_id, cos_dec_id, dec_val=None):
        if self.wc_fisher:
            cov_ra_cos_dec = self.cov[ra_id,cos_dec_id]
            if dec_val is None:
                self.errs['sky_area_90'] = edh.sky_area_90(self.errs['ra'],self.errs['cos_dec'],cov_ra_cos_dec,dec_val)
            else:
                self.errs['sky_area_90'] = edh.sky_area_90(self.errs['ra'],self.errs['dec'],cov_ra_cos_dec,dec_val)


    ###
    #-----IO methods-----
    def print_detector(self,print_format=1):
        if print_format:
            sepl='-----------------------------------------------------------------------------------'
            print()
            print(sepl)
            print('Printing detector.')
            print(sepl)
            print()
        for key,value in vars(self).items():
            if type(value) == dict:
                print('Key: ',key)
                for key in value.keys():
                    print('',key)
                    print('',value[key])
                print()
            elif value is not None:
                print('Key: ',key)
                print(value)
                print()
        if print_format:
            print(sepl)
            print('Printing detector done.')
            print(sepl)
            print()


###
#-----Convert the evaluated derivatives according to conv_cos, conv_log-----
def get_conv_inj_params_deriv_variables(c_quants, inj_params, deriv_variables):
    if c_quants == {}:
        return inj_params, deriv_variables

    else:
        for o_key in c_quants:
            c_key = c_quants[o_key][0]
            c_val = c_quants[o_key][1]

            if c_key not in deriv_variables:
                deriv_variables[deriv_variables.index(o_key)] = c_key
            inj_params[c_key] = c_val

        return inj_params, deriv_variables

def get_conv_del_eval_dic(del_eval_dic, params_dic, conv_cos, conv_log, deriv_symbs_string):
    if conv_cos is None and conv_log is None:
        return del_eval_dic, {}
    else:
        conv_dic = {}
        c_quants = {}

        deriv_variables = deriv_symbs_string.split(' ')
        for deriv in del_eval_dic:
            for deriv_var in deriv_variables:
                if deriv_var in deriv:
                    key = deriv_var
                    break

            c_key = None
            if conv_cos is not None and key in conv_cos:
                c_deriv, c_key, c_val = edh.convert_to_cos_derivative(del_eval_dic[deriv],key,params_dic[key])
            elif conv_log is not None and key in conv_log:
                c_deriv, c_key, c_val = edh.convert_to_log_derivative(del_eval_dic[deriv],key,params_dic[key])
            else:
                conv_dic[deriv] = del_eval_dic[deriv]

            if c_key is not None:
                c_quants[key] = (c_key, c_val)

                prefix,suffix = deriv.split(key)
                n_key = prefix + c_key + suffix
                conv_dic[n_key] = c_deriv

        return conv_dic, c_quants
