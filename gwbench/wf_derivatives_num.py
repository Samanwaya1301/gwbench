'''This module contains two methods that calculate numerical derivatives.
'''

import numdifftools as nd

import gwbench.basic_functions as bfs
import gwbench.wf_manipulations as wfm

def part_deriv_hf_func(hf, symbols_list, deriv_symbs_list, f, params_dic, pl_cr=0, compl=1, label='hf', step=1e-7, method='central', order=2, n=1):

    if 'f' in symbols_list:
        symbols_list.remove('f')

    params_list = []
    for param in symbols_list:
        params_list.append(params_dic[param])

    deriv_params_list = []
    for param in deriv_symbs_list:
        deriv_params_list.append(params_dic[param])

    def hf_of_deriv_params(f, *deriv_params_list):
        tmp_list = []
        for i,el in enumerate(symbols_list):
            if el in deriv_symbs_list:
                tmp_list.append(deriv_params_list[deriv_symbs_list.index(el)])
            else:
                tmp_list.append(params_list[i])

        return hf(f,*tmp_list)

    del_hf = part_deriv(hf_of_deriv_params, f, deriv_params_list, pl_cr, compl, step, method, order, n)
    if pl_cr:
        del_hf_dic = {}

        if len(deriv_symbs_list) == 1:
            key_string = 'del_'+deriv_symbs_list[0]+'_'+label+'p'
            del_hf_dic[key_string] = del_hf[0]

            key_string = 'del_'+deriv_symbs_list[0]+'_'+label+'c'
            del_hf_dic[key_string] = del_hf[1]

        else:
            for i,name in enumerate(deriv_symbs_list):
                key_string = 'del_'+name+'_'+label+'p'
                del_hf_dic[key_string] = del_hf[0][:,i]

                key_string = 'del_'+name+'_'+label+'c'
                del_hf_dic[key_string] = del_hf[1][:,i]

    else:
        del_hf_dic = {}

        if len(deriv_symbs_list) == 1:
            key_string = 'del_'+deriv_symbs_list[0]+'_'+label
            del_hf_dic[key_string] = del_hf

        else:
            for i,name in enumerate(deriv_symbs_list):
                key_string = 'del_'+name+'_'+label
                del_hf_dic[key_string] = del_hf[:,i]

    return del_hf_dic


def part_deriv(func, f, params_list, pl_cr=0, compl=None, step=1e-7, method='central', order=2, n=1):
    if pl_cr:
        if compl:
            def amp_pha_of_hfpc_compl(f, *params_list):
                return wfm.pl_cr_to_amp_pha(*func(f, *params_list))

            amp_pl, pha_pl, amp_cr, pha_cr = amp_pha_of_hfpc_compl(f, *params_list)

            del_amp_pl = part_deriv_ndGradient(amp_pha_of_hfpc_compl, f, params_list, 0, step, method, order, n)
            del_pha_pl = part_deriv_ndGradient(amp_pha_of_hfpc_compl, f, params_list, 1, step, method, order, n)
            del_amp_cr = part_deriv_ndGradient(amp_pha_of_hfpc_compl, f, params_list, 2, step, method, order, n)
            del_pha_cr = part_deriv_ndGradient(amp_pha_of_hfpc_compl, f, params_list, 3, step, method, order, n)

            del_hfp = wfm.z_deriv_from_amp_pha(amp_pl, pha_pl, del_amp_pl, del_pha_pl)
            del_hfc = wfm.z_deriv_from_amp_pha(amp_cr, pha_cr, del_amp_cr, del_pha_cr)

            return del_hfp, del_hfc

        else:
            def amp_pha_of_hfpc_real(f, *params_list):
                return wfm.amp_pha_from_re_im(*func(f, *params_list))

            amp, pha = amp_pha_of_hfpc_real(f, *params_list)
            del_amp = part_deriv_ndGradient(amp_pha_of_hfpc_real, f, params_list, 0, step, method, order, n)
            del_pha = part_deriv_ndGradient(amp_pha_of_hfpc_real, f, params_list, 1, step, method, order, n)

            return wfm.re_im_from_z(wfm.z_deriv_from_amp_pha(amp, pha, del_amp, del_pha))

    else:
        if compl:
            def amp_pha_of_hf(f, *params_list):
                return wfm.amp_pha_from_z(func(f, *params_list))

            amp, pha = amp_pha_of_hf(f, *params_list)
            del_amp = part_deriv_ndGradient(amp_pha_of_hf, f, params_list, 0, step, method, order, n)
            del_pha = part_deriv_ndGradient(amp_pha_of_hf, f, params_list, 1, step, method, order, n)

            return wfm.z_deriv_from_amp_pha(amp, pha, del_amp, del_pha)

        else:
            return part_deriv_ndGradient(func, f, params_list, None, step, method, order, n)


def part_deriv_ndGradient(func, f, params_list=None, funcid=None, step=1e-7, method='central', order=2, n=1):
    def wraps(x):
        if funcid == None:
            return func(f,*x)
        else:
            return func(f,*x)[funcid]

    if params_list == None:
        return nd.Gradient(wraps, step=step, method=method, order=order, n=n)
    else:
        return nd.Gradient(wraps, step=step, method=method, order=order, n=n)(params_list)
