import os
import sys

import dill

import gwbench.antenna_pattern_np as ant_pat_np
import gwbench.antenna_pattern_sp as ant_pat_sp
import gwbench.basic_functions as bfs
import gwbench.wf_derivatives_num as wfd_num
import gwbench.wf_derivatives_sym as wfd_sym

lambdified_functions_path = os.path.join(os.getcwd(),'lambdified_functions')
ant_pat_symbs_string = 'f Mc tc ra dec psi gmst0'


def calc_det_responses_num(loc, wf, deriv_symbs_string, f_arr, params_dic, use_rot=1, label='hf', step=1e-7, method='central', order=2, n=1):

    wf_symbs_list = wf.wf_symbs_string.split(' ')
    deriv_symbs_list = deriv_symbs_string.split(' ')

    if 'f' in wf_symbs_list:
        wf_symbs_list.remove('f')
    if 'f' in deriv_symbs_list:
        deriv_symbs_list.remove('f')

    if loc == None:
        wf_params_list = list(bfs.get_sub_dict(params_dic,wf_symbs_list).values())

        def pc_func(f_arr,*wf_params_list):
            wf_list = []
            for i,el in enumerate(wf_symbs_list):
                wf_list.append(wf_params_list[wf_symbs_list.index(el)])
            return wf.eval_np_func(f_arr, wf_list)

        return wfd_num.part_deriv_hf_func(pc_func, wf_symbs_list, deriv_symbs_list, f_arr, params_dic, pl_cr=1, compl=1, label=label, step=step, method=method, order=order, n=n)

    else:
        ap_symbs_list = ant_pat_symbs_string.split(' ')
        if 'f' in ap_symbs_list:
            ap_symbs_list.remove('f')

        dr_symbs_list = bfs.reduce_symbols_strings(wf.wf_symbs_string,ant_pat_symbs_string).split(' ')
        dr_params_list = list(bfs.get_sub_dict(params_dic,dr_symbs_list).values())

        def dr_func(f_arr,*dr_params_list):
            wf_list = []
            for i,el in enumerate(wf_symbs_list):
                wf_list.append(dr_params_list[dr_symbs_list.index(el)])

            ap_list = []
            for i,el in enumerate(ap_symbs_list):
                ap_list.append(dr_params_list[dr_symbs_list.index(el)])

            hfp, hfc = wf.eval_np_func(f_arr, wf_list)
            Fp, Fc, Flp = ant_pat_np.antenna_pattern_and_loc_phase_fac(f_arr,*ap_list,loc,use_rot)

            return Flp * (hfp * Fp + hfc * Fc)

        return wfd_num.part_deriv_hf_func(dr_func, dr_symbs_list, deriv_symbs_list, f_arr, params_dic, pl_cr=0, compl=1, label=label, step=step, method=method, order=order, n=n)



def generate_det_responses_sym(wf,deriv_symbs_string,locs=None,use_rot=1):

    hfpc = wf.get_sp_expr()

    responses = {}
    responses['pl_cr'] = hfpc

    if locs is None:
        locs = ant_pat_np.locs
    else:
        for loc in locs:
            if loc not in ant_pat_np.locs:
                exit_str =('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'+
                          f'Could not find the lambdified function file: {loc}\n'+
                           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                sys.exit(exit_str) 
                
    for loc in locs:
        responses[loc] = ant_pat_sp.detector_response_expr(hfpc[0],hfpc[1],loc,use_rot)

    for key in responses.keys():
        if key == 'pl_cr':
            print('Calculating the derivatives of the plus/cross polarizations.')
            wf_deriv_symbs_string = bfs.remove_symbols(deriv_symbs_string,wf.wf_symbs_string)
            deriv_dic = wfd_sym.part_deriv_hf_expr(responses[key],wf.wf_symbs_string,wf_deriv_symbs_string,pl_cr=1)
            deriv_dic['variables'] = wf.wf_symbs_string
            deriv_dic['deriv_variables'] = wf_deriv_symbs_string

            file_name = 'par_deriv_WFM_'+wf.wf_model_name+'_VAR_'+wf_deriv_symbs_string.replace(' ', '_')+'_DET_'+key+'.dat'

        else:
            print('Calculating the derivatives of the detector response for detector: ' + key)
            symbols_string = bfs.reduce_symbols_strings(wf.wf_symbs_string,ant_pat_symbs_string)
            deriv_dic = wfd_sym.part_deriv_hf_expr(responses[key],symbols_string,deriv_symbs_string)
            deriv_dic['variables'] = symbols_string
            deriv_dic['deriv_variables'] = deriv_symbs_string

            file_name = 'par_deriv_WFM_'+wf.wf_model_name+'_VAR_'+deriv_symbs_string.replace(' ', '_')+'_DET_'+key+'.dat'

        file_name = os.path.join(lambdified_functions_path,file_name)

        if not os.path.exists(lambdified_functions_path):
            os.makedirs(lambdified_functions_path)

        with open(file_name, "wb") as fi:
            dill.dump(deriv_dic, fi, recurse=True)

    print('Done.')
    return


def load_det_response_sym(det_name, wf_model_name, deriv_symbs_string, return_bin=0):
    file_name = 'par_deriv_WFM_'+wf_model_name+'_VAR_'+deriv_symbs_string.replace(' ', '_')+'_DET_'+det_name+'.dat'
    file_name = os.path.join(lambdified_functions_path,file_name)

    try:
        with open(file_name, "rb") as fi:
            if return_bin:
                return fi.read()
            else:
                return dill.load(fi)
    except FileNotFoundError:
        exit_str = ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n' + 
                   f'Could not find the lambdified function file: {file_name}\n' +
                    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        sys.exit(exit_str) 
