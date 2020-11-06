# Copyright (C) 2020  Ssohrab Borhanian
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


"""This module handles the benchmarking of graviational waveforms observed by a network of detectors.

"""

from copy import copy
from multiprocessing import Pool

import dill
import numpy as np

import gwbench.antenna_pattern_np as ant_pat_np
import gwbench.basic_functions as bfs
import gwbench.detector_class as dc
import gwbench.detector_responses as dr
import gwbench.err_deriv_handling as edh
import gwbench.fisher_analysis_tools as fat
import gwbench.wf_class as wfc

class Network:

    ###
    #-----Init methods-----
    def __init__(self,network_spec=None):
        ##---initialize network object-----
        if network_spec is None:
            #-----network and detectors
            # network label
            self.label = None
            # detector labels list
            self.det_keys = None
            # list of detectors in network
            self.detectors = None

        elif isinstance(network_spec, str):
            #-----network and detectors
            self.set_network_and_detectors_from_label(network_spec)
        elif isinstance(network_spec, list) or isinstance(network_spec, tuple):
            #-----network and detectors
            self.set_network_and_detectors_from_key_list(network_spec)

        #-----injection and waveform quantities-----
        # frequency array
        self.f = None
        # dictionary of injection parameters
        self.inj_params = None
        # derivative variables - symbs_string and list
        self.deriv_symbs_string = None
        self.deriv_variables = None
        # waveform
        self.wf = None
        # antenna pattern and location phase factor symbs_string
        self.ap_symbs_string = ant_pat_np.ap_symbs_string

        #-----analysis settings-----
        # list of inj_params to convert to cos, ln versions
        self.conv_cos = None
        self.conv_log = None
        # use f-dependent gmst (SPA) in antenna patterns
        self.use_rot = None

        #-----waveform polarizations-----
        # plus/cross polarizations
        self.hfp = None
        self.hfc = None
        # derivative dictionary for polarizations
        self.del_hfpc = None
        # sympy expressions of derivative dictionary for polarizations
        self.del_hfpc_expr = None

        #-----network SNR-----
        # SNR, SNR^2 calculated from hf
        self.snr = None
        self.snr_sq = None

        #-----network errors-----
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
    #
    # it is best practice to always change the instance variables using these setter methods
    #
    def set_network_and_detectors_from_label(self, network_label):
        #-----network and detectors
        self.label = network_label
        self.det_keys = read_det_keys_from_label(network_label)
        self.detectors = []
        for det_key in self.det_keys:
            self.detectors.append(dc.Detector(det_key))

    def set_network_and_detectors_from_key_list(self, det_key_list):
        if isinstance(det_key_list, tuple):
            tmp_list = []
            for tec,loc in zip(det_key_list[0],det_key_list[1]):
                tmp_list.append(tec + '_' + loc)
            det_key_list = tmp_list

        #-----network and detectors
        self.label = '..'.join(det_key_list)
        self.det_keys = det_key_list
        self.detectors = []
        for det_key in self.det_keys:
            self.detectors.append(dc.Detector(det_key))

    def set_wf_vars(self, wf_model_name, wf_other_var_dic=None):
        self.wf = wfc.Waveform(wf_model_name, wf_other_var_dic)

    def set_net_vars(self, f=None, inj_params=None, deriv_symbs_string=None, conv_cos=None, conv_log=None, use_rot=None):
        if f is not None:
            self.f = copy(f)
            if self.detectors is not None:
                for det in self.detectors:
                    det.set_f(self.f)
        if inj_params is not None:
            self.inj_params = copy(inj_params)
        if deriv_symbs_string is not None:
            self.deriv_symbs_string = copy(deriv_symbs_string)
            self.deriv_variables = deriv_symbs_string.split(' ')
        if conv_cos is not None:
            self.conv_cos = copy(conv_cos)
        if conv_log is not None:
            self.conv_log = copy(conv_log)
        if use_rot is not None:
            self.use_rot = copy(use_rot)

    ##
    #-----Resetter methods for instance variables-----
    def reset_ant_pat_lpf_psds(self):
        for det in self.detectors:
            det.Fp = None
            det.Fc = None
            det.Flp = None
            det.psd = None

    def reset_wf_polarizations(self):
        self.hfp = None
        self.hfc = None
        self.del_hfpc = None
        self.del_hfpc_expr = None

    def reset_det_responses(self):
        for det in self.detectors:
            det.hf = None
            det.del_hf = None
            det.del_hf_expr = None

    def reset_snrs(self):
        self.snr = None
        self.snr_sq = None
        for det in self.detectors:
            det.snr = None
            det.snr_sq = None
            det.d_snr_sq = None

    def reset_errors(self):
        self.fisher = None
        self.cond_num = None
        self.wc_fisher = None
        self.cov = None
        self.inv_err = None
        self.errs = None
        for det in self.detectors:
            det.fisher = None
            det.cond_num = None
            det.wc_fisher = None
            det.cov = None
            det.inv_err = None
            det.errs = None


    ###
    #-----Getters-----
    def get_detector(self,det_key):
        return self.detectors[self.det_keys.index(det_key)]

    def get_snrs_errs_cov_fisher_inv_err_for_key(self,key='network'):
        if key == 'network':
            out_obj = self
        else:
            out_obj = self.detectors[self.det_keys.index(key)]
        snr = out_obj.snr
        return snr, out_obj.errs, out_obj.cov, out_obj.fisher, out_obj.inv_err


    ###
    #-----PSDs and antenna patterns-----
    def setup_ant_pat_lpf_psds(self, F_lo=-np.inf, F_hi=np.inf, psd_file_dict=None):
        for det in self.detectors:
            det.setup_ant_pat_lpf_psds(self.inj_params, self.use_rot, F_lo, F_hi, psd_file_dict)
        print('Antenna patterns, LPFs, and PSDs loaded.')

    def setup_psds(self, F_lo=-np.inf, F_hi=np.inf, psd_file_dict=None):
        for det in self.detectors:
            det.setup_psds(F_lo, F_hi, psd_file_dict)
        print('PSDs loaded.')

    def setup_ant_pat_lpf(self):
        for det in self.detectors:
            det.setup_ant_pat_lpf(self.inj_params, self.use_rot)
        print('Antenna patterns and LPFs loaded.')


    ###
    #-----Waveform polarizations-----
    def calc_wf_polarizations(self):
        self.hfp, self.hfc = self.wf.eval_np_func(self.f,bfs.get_sub_dict(self.inj_params,self.wf.wf_symbs_string))
        print('Polarizations calculated.')

    def calc_wf_polarizations_derivs_num(self, step=1e-7, method='central', order=2, n=1):
        print('Calculate numeric derivatives of polarizations.')
        self.calc_wf_polarizations()
        wf_deriv_symbs_string = bfs.remove_symbols(self.deriv_symbs_string,self.wf.wf_symbs_string)
        self.del_hfpc = dr.calc_det_responses_num(None,self.wf,wf_deriv_symbs_string,self.f,self.inj_params,self.use_rot,'hf',step,method,order,n)
        self.del_hfpc, c_quants = dc.get_conv_del_eval_dic(self.del_hfpc, self.inj_params, self.conv_cos, self.conv_log, self.deriv_symbs_string)
        self.inj_params, self.deriv_variables = dc.get_conv_inj_params_deriv_variables(c_quants, self.inj_params, self.deriv_variables)
        print('Numeric derivatives of polarizations calculated.')

    def load_wf_polarizations_derivs_sym(self, return_bin=0):
        wf_deriv_symbs_string = bfs.remove_symbols(self.deriv_symbs_string,self.wf.wf_symbs_string)
        self.del_hfpc_expr = dr.load_det_response_sym('pl_cr', self.wf.wf_model_name, wf_deriv_symbs_string, return_bin)
        print('Lambdified polarizations loaded.')

    def calc_wf_polarizations_derivs_sym(self):
        print('Evaluate polarizations.')
        self.calc_wf_polarizations()
        self.del_hfpc = {}
        for deriv in self.del_hfpc_expr:
            if deriv in ('variables','deriv_variables'): continue
            self.del_hfpc[deriv] = self.del_hfpc_expr[deriv](self.f, **bfs.get_sub_dict(self.inj_params, self.del_hfpc_expr['variables']))
        self.del_hfpc, c_quants = dc.get_conv_del_eval_dic(self.del_hfpc, self.inj_params, self.conv_cos, self.conv_log, self.deriv_symbs_string)
        self.inj_params, self.deriv_variables = dc.get_conv_inj_params_deriv_variables(c_quants, self.inj_params, self.deriv_variables)
        print('Lambdified polarizations evaluated.')


    ###
    #-----Detector responses-----
    def calc_det_responses(self):
        for det in self.detectors:
            det.calc_det_responses(self.wf,self.inj_params)
        print('Detector responses calculated.')

    def calc_det_responses_derivs_num(self, step=1e-7, method='central', order=2, n=1):
        print('Calculate numeric derivatives of detector responses.')
        for det in self.detectors:
            det.calc_det_responses_derivs_num(self.inj_params, self.deriv_variables, self.wf, self.deriv_symbs_string, self.conv_cos, self.conv_log, self.deriv_symbs_string, self.use_rot, step, method, order, n)
        print('Numeric derivatives of detector responses calculated.')

    def load_det_responses_derivs_sym(self, return_bin=0):
        for det in self.detectors:
            det.load_det_responses_derivs_sym(self.wf.wf_model_name, self.deriv_symbs_string, return_bin)
        print('Lambdified detector responses loaded.')

    def calc_det_responses_derivs_sym(self):
        print('Evaluate lambdified detector responses.')
        for det in self.detectors:
            det.calc_det_responses_derivs_sym(self.wf, self.inj_params, self.deriv_variables, self.conv_cos, self.conv_log, self.deriv_symbs_string)
        print('Lambdified detector responses evaluated.')


    ###
    #-----SNR calculations-----
    def calc_snrs_det_responses(self, only_net=0, df=None):
        print('Calculate SNRs.')
        self.snr_sq = 0
        for det in self.detectors:
             self.snr_sq += det.calc_snrs_det_responses(only_net, df)
        self.snr = np.sqrt(self.snr_sq)
        print('SNRs calculated.')

    def calc_snr_sq_integrand_det_responses(self):
        print('Calculate SNR integrands.')
        for det in self.detectors:
            det.calc_snr_sq_integrand_det_responses()
        print('SNR integrands calculated.')


    ###
    #-----Error calculation and Fisher analysis-----
    def calc_errors(self, cond_sup=1e15, by_element=0, only_net=0, df=None):
        print('Calculate errors (Fisher & cov matrices).')
        #-----calculate the error matrices: Fisher and Cov-----
        self.fisher = 0
        for det in self.detectors:
            self.fisher += det.calc_error_mats(only_net, df, cond_sup)
        self.wc_fisher, self.cond_num = fat.check_well_conditioned(self.fisher,cond_sup)
        self.cov = fat.calc_cov_from_fisher(self.fisher,self.wc_fisher)
        #-----calculate the max inversion errors of the detectors and network-----
        self.inv_err = fat.inv_err_from_fisher_cov(self.fisher,self.cov,by_element)
        if not only_net:
            for det in self.detectors:
                det.calc_inv_err(by_element)
        #-----calculate the absolute errors of the various variables-----
        self.errs = fat.get_errs_from_cov(self.cov,self.deriv_variables)
        if not only_net:
            for det in self.detectors:
                det.calc_errs(self.deriv_variables)
        print('Errors calculated.')

    def calc_sky_area_90(self,only_net=0):
        print('Calculate 90% sky area.')
        if 'ra' in self.deriv_variables and ('cos_dec' in self.deriv_variables or 'dec' in self.deriv_variables):
            ra_id = self.deriv_variables.index('ra')
            if 'cos_dec' in self.deriv_variables:
                cos_dec_id  = self.deriv_variables.index('cos_dec')
                dec_val     = None
            else:
                cos_dec_id  = self.deriv_variables.index('dec')
                dec_val     = self.inj_params['dec']
            if self.wc_fisher:
                cov_ra_cos_dec = self.cov[ra_id,cos_dec_id]
                if dec_val is None:
                    self.errs['sky_area_90'] = edh.sky_area_90(self.errs['ra'],self.errs['cos_dec'],cov_ra_cos_dec,dec_val)
                else:
                    self.errs['sky_area_90'] = edh.sky_area_90(self.errs['ra'],self.errs['dec'],cov_ra_cos_dec,dec_val)
            if not only_net:
                for det in self.detectors:
                    det.calc_sky_area_90(ra_id,cos_dec_id,dec_val)
            print('Sky area calculated.')
        else:
            print('Nothing done due to missing of either RA or COS_DEC (DEC) errors.')


    ###
    #-----IO methods-----
    def save_network(self,filename_path):
        '''Save the network under the given path using *dill*.'''
        with open(filename_path, "wb") as fi:
            dill.dump(self, fi, recurse=True)
        print('Network pickled.')
        return

    def load_network(self,filename_path):
        '''Loading the network from the given path using *dill*.'''
        with open(filename_path, "rb") as fi:
            self = dill.load(fi)
        print('Network loaded.')
        return network

    def print_network(self):
        sepl='--------------------------------------------------------------------------------------'
        print()
        print(sepl)
        print('Printing network.')
        print(sepl)
        print()
        for key,value in vars(self).items():
            if type(value) == dict:
                print('Key: ',key)
                for kkey in value.keys():
                    print('',kkey)
                    print('',value[kkey])
                print()
            elif value is not None:
                if key == 'wf':
                    print('Key: ',key)
                    for kkey,vvalue in vars(value).items():
                        print('',kkey.ljust(16,' '),'  ',vvalue)
                    print()
                else:
                    print('Key: ',key)
                    print(value)
                    print()
        print(sepl)
        print('Printing network done.')
        print(sepl)
        print()

    def print_detectors(self):
        sepl='--------------------------------------------------------------------------------------'
        sepl1='-------------------------------------------'
        print()
        print(sepl)
        print('Printing detectors.')
        print(sepl)
        for det in self.detectors:
            print(sepl1)
            print(det.det_key)
            print(sepl1)
            det.print_detector(0)
        print(sepl)
        print('Printing detectors done.')
        print(sepl)
        print()

    ###
    #-----Dealing with several networks-----
    def get_det_responses_psds_from_locs_tecs(self,loc_net,tec_net,sym_derivs=0):

        self.inj_params = loc_net.inj_params
        self.deriv_variables = loc_net.deriv_variables
        self.f = loc_net.f

        for i,det in enumerate(self.detectors):
            tec_det = tec_net.get_detector(det.tec+'_loc')
            det.f = tec_det.f
            det.psd = tec_det.psd

            f_lo = det.f[0]
            f_hi = det.f[-1]
            ids = np.logical_and(self.f>=f_lo,self.f<=f_hi)

            loc_det = loc_net.get_detector('tec_'+det.loc)
            if sym_derivs:
                det.del_hf_expr = copy(loc_det.del_hf_expr)
            det.hf = copy(loc_det.hf[ids])
            det.del_hf = copy(loc_det.del_hf)
            det.deriv_variables = copy(loc_net.deriv_variables)
            det.inj_params = copy(loc_net.inj_params)

            for deriv in det.del_hf:
                det.del_hf[deriv] = det.del_hf[deriv][ids]

        print('Detector responses transferred.')


###
#-----Dealing with several networks-----
def unique_tecs(network_labels,f,F_lo=-np.inf,F_hi=np.inf,psd_file_dict=None):
    print('Calculate PSDs for unique detector technologies.')

    # initialize empty network
    tec_net = Network()
    # get the detector keys
    tec_net.det_keys = []

    # find unique technologies
    for network_label in network_labels:
        network = Network(network_label)
        for det in network.detectors:
            tec_net.det_keys.append(det.tec)
    tec_net.det_keys = list(dict.fromkeys(tec_net.det_keys))

    # make them into fake detector keys
    for i,tec in enumerate(tec_net.det_keys):
        tec_net.det_keys[i] = tec+'_loc'
    # initialize fake detectors
    tec_net.detectors = []
    for det_key in tec_net.det_keys:
        tec_net.detectors.append(dc.Detector(det_key))

    # get PSDs
    tec_net.set_net_vars(f=f)
    tec_net.setup_psds(F_lo,F_hi,psd_file_dict)

    print('PSDs for unique detector technologies calculated.')
    return tec_net


def unique_locs_det_responses(network_labels,f,inj_params,deriv_symbs_string,wf_model_name,wf_other_var_dic=None,conv_cos=None,conv_log=None,use_rot=1,num_cores=None,step=None,method=None,order=None,n=None):
    print('Evaluate lambdified detector responses for unique locations.')

    # initialize empty network
    loc_net = Network()
    # get the detector keys
    loc_net.det_keys = []

    # find unique locations
    for network_label in network_labels:
        network = Network(network_label)
        for det in network.detectors:
            loc_net.det_keys.append(det.loc)
    loc_net.det_keys = list(dict.fromkeys(loc_net.det_keys))

    # make them into fake detectpr keys
    for i,loc in enumerate(loc_net.det_keys):
        loc_net.det_keys[i] = 'tec_'+loc
    # initialize fake detectors
    loc_net.detectors = []
    for det_key in loc_net.det_keys:
        loc_net.detectors.append(dc.Detector(det_key))
    # set all the other necessary variables
    loc_net.set_net_vars(f=f, inj_params=inj_params, deriv_symbs_string=deriv_symbs_string, conv_cos=conv_cos, conv_log=conv_log, use_rot=use_rot)
    loc_net.set_wf_vars(wf_model_name,wf_other_var_dic)

    # setup Fp, Fc, and Flp and calculate the detector responses
    loc_net.setup_ant_pat_lpf()
    loc_net.calc_det_responses()

    if step is None:
        print('Loading the lamdified functions.')
        loc_net.load_det_responses_derivs_sym(return_bin = 1)
        print('Loading done.')


    print('Starting evaluation.')
    if num_cores is None:
        if step is None:
            for det in loc_net.detectors:
                det.del_hf, c_quants = eval_loc_sym(det.loc,det.del_hf_expr,deriv_symbs_string,f,inj_params,conv_cos,conv_log)
        else:
            for det in loc_net.detectors:
                det.del_hf, c_quants = eval_loc_num(det.loc,loc_net.wf,deriv_symbs_string,f,inj_params,conv_cos,conv_log,use_rot,step,method,order,n)

        loc_net.inj_params, loc_net.deriv_variables = dc.get_conv_inj_params_deriv_variables(c_quants, loc_net.inj_params, loc_net.deriv_variables)

    else:
        pool = Pool(num_cores)
        if step is None:
            arg_tuple_list = [(det.loc,det.del_hf_expr,deriv_symbs_string,f,inj_params,conv_cos,conv_log) for det in loc_net.detectors]
            result = pool.starmap_async(eval_loc_sym, arg_tuple_list)
            result.wait()
        else:
            arg_tuple_list = [(det.loc,loc_net.wf,deriv_symbs_string,f,inj_params,conv_cos,conv_log,use_rot,step,method,order,n) for det in loc_net.detectors]
            result = pool.starmap_async(eval_loc_num, arg_tuple_list)
            result.wait()

        for det, (del_hf,c_quants) in zip(loc_net.detectors, result.get()):
            det.del_hf = del_hf

        loc_net.inj_params, loc_net.deriv_variables = dc.get_conv_inj_params_deriv_variables(c_quants, loc_net.inj_params, loc_net.deriv_variables)

    if step is None:
        for det in loc_net.detectors:
            det.del_hf_expr = dill.loads(det.del_hf_expr)

    print('Lambdified detector responses for unique locations evaluated.')
    return loc_net

def eval_loc_sym(loc,del_hf_expr,deriv_symbs_string,f,inj_params,conv_cos,conv_log):
    print(' ',loc)
    del_hf = {}
    del_hf_expr = dill.loads(del_hf_expr)
    for deriv in del_hf_expr:
        if deriv in ('variables','deriv_variables'): continue
        del_hf[deriv] = del_hf_expr[deriv](f,**bfs.get_sub_dict(inj_params,del_hf_expr['variables']))
    return dc.get_conv_del_eval_dic(del_hf,inj_params,conv_cos,conv_log, deriv_symbs_string)

def eval_loc_num(loc,wf,deriv_symbs_string,f,inj_params,conv_cos,conv_log,use_rot,step,method,order,n):
    print(' ',loc)
    del_hf = dr.calc_det_responses_num(loc,wf,deriv_symbs_string,f,inj_params,use_rot,'hf',step,method,order,n)
    return dc.get_conv_del_eval_dic(del_hf,inj_params,conv_cos,conv_log, deriv_symbs_string)


###
#-----Get detector list for strings-----
def read_det_keys_from_label(network_label):
    det_keys = []

    ##-----read the network label and find all detectors-----
    keys = list(network_label)

    # - in network list means that all 2G detectors up to that index are to be
    # taken at aLIGO sensitivity
    aLIGO = int('-' in keys)
    if aLIGO: aLIGO_id = keys.index('-')

    # + in network list means that all 2G detectors up to that index are to be
    # taken at A+ sensitivity
    a_pl = int('+' in keys)
    if a_pl: a_pl_id = keys.index('+')

    # v in network list means that all 2G detectors up to that index are to be
    # taken at Voyager sensitivity
    voy = int('v' in keys)
    if voy:
        voy_id = keys.index('v')
        tmp = int(keys[voy_id+1] == 'p')
        voy_pmo = tmp * 'PMO' + (1-tmp) * 'CBO'

    # find out which locations with which PSDs are in the network
    for loc in ant_pat_np.locs:
        if loc in keys:
            loc_id = keys.index(loc)

            if ant_pat_np.check_loc_gen(loc) == '2G':
                if aLIGO and loc_id < aLIGO_id:
                    name = 'aLIGO_'+loc
                elif a_pl and loc_id < a_pl_id:
                    if loc == 'V':
                        name = 'V+_'+loc
                    elif loc == 'K':
                        name = 'K+_'+loc
                    else:
                        name = 'A+_'+loc
                elif voy and loc_id < voy_id:
                    name = 'Voyager-{}_{}'.format(voy_pmo,loc)

            else:
                ce_a = int(keys[loc_id+1] == 'a') # 0 for i, 1 for a - CE1 as i, CE2 as a
                ce_arm = int(keys[loc_id+2])*10  # arm length (n for n*10km)
                tmp = int(keys[loc_id+3] == 'p')
                ce_pmo = tmp * 'PMO' + (1-tmp) * 'CBO'
                name = 'CE%i-%i-{}_{}'.format(ce_pmo,loc)%(ce_a+1,ce_arm)

            det_keys.append(name)

    # add 3 ET detectors
    if 'E' in keys:
        for name in ['ET_ET1','ET_ET2','ET_ET3']:
            det_keys.append(name)

    return det_keys
