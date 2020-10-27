import numpy as np

import gwbench.snr as snr_mod

def fisher_cov_matrix(del_vs_f_list,psd,f,only_fisher=0,df=None,cond_sup=1e15):
    n = len(del_vs_f_list)
    fisher = np.zeros((n,n))

    for i in np.arange(n):
        fisher[i,i] = snr_mod.scalar_product_freq_array(del_vs_f_list[i],del_vs_f_list[i],psd,f,df)
        for j in np.arange(i+1,n):
           fisher[i,j] = snr_mod.scalar_product_freq_array(del_vs_f_list[i],del_vs_f_list[j],psd,f,df)
           fisher[j,i] = fisher[i,j]

    wc_fisher, cond_num = check_well_conditioned(fisher,cond_sup)
    # return cov=None, if Fisher not well conditioned OR if we only fisher is wanted
    cov = calc_cov_from_fisher(fisher,(wc_fisher and not only_fisher))
    return fisher, cov, wc_fisher, cond_num

def calc_cond_number(fisher):
    EWs,_ = np.linalg.eig(fisher)
    return np.amax(np.abs(EWs))/np.amin(np.abs(EWs))

def check_well_conditioned(fisher,cond_sup=1e15):
    if cond_sup is None: cond_sup = np.inf
    cond_num = calc_cond_number(fisher)
    return ( cond_num < cond_sup), cond_num

def calc_cov_from_fisher(fisher,wc_fisher):
    if wc_fisher: return np.linalg.inv(fisher)
    else:         return None

def inv_err_from_fisher_cov(fisher,cov,by_element=0):
    if cov is None:
        return None
    else:
        ident = np.identity(cov.shape[0])
        res = np.abs(np.matmul(fisher,cov)-ident)
        if by_element: return np.array([np.maximum(np.amax(res[:,i]),np.amax(res[i])) for i in range(cov.shape[0])])
        else:          return np.amax(res)

def get_errs_from_cov(cov,deriv_variables):
    if cov is None:
        return None
    else:
        errs = {}
        for i,name in enumerate(deriv_variables):
            errs[name] = np.sqrt(np.abs(cov[i,i]))
        return errs
