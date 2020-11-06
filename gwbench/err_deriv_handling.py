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


import numpy as np
from scipy.integrate import simps
from scipy.optimize import toms748

PI = np.pi

#-----derivative and error manipulations-----
def convert_to_cos_derivative(deriv, param_key, param_val):
    c_deriv = (-1./np.sin(param_val)) * deriv
    c_param_key = 'cos_'+param_key
    c_param_val = np.cos(param_val)
    return c_deriv, c_param_key, c_param_val

def convert_to_log_derivative(deriv, param_key, param_val):
    c_deriv = param_val * deriv
    c_param_key = 'log_'+param_key
    c_param_val = np.log(param_val)
    return c_deriv, c_param_key, c_param_val

def dim_err_to_rel_err(err, param_val, param_key=None, param_kind=None):
    if param_key in ('M','Mc','Dl','DL','tc') or param_kind == 'dim':
        return np.abs(err/param_val)
    else:
        return err

def one_sigma_to_percent_error(percent, sigma):
    sigma_orders = [0, 1, 2, 3, 4, 5, 6]
    sigma_percents = [0., 68.2689492137, 95.4499736104, 99.7300203937, 99.9936657516, 99.9999426697, 99.9999998027]
    for i,sigma_order in enumerate(sigma_orders):
        if percent < sigma_percents[i]: break

    def func(error,percent,sigma,sigma_order):
        x = np.linspace(0,error,sigma_order*100)
        return percent/100/2 - simps(np.exp(-x**2/2/sigma**2)/np.sqrt(2*np.pi)/sigma,x)

    if sigma_order == 1:
        return toms748(func,1e-2*sigma,sigma_order*sigma,args=(percent,sigma,sigma_order))
    else:
        return toms748(func,(sigma_order-1)*sigma,sigma_order*sigma,args=(percent,sigma,sigma_order))

#-----sky area calculations-----
def sky_area_90(ra_err, cos_dec_err, cov_ra_cos_dec, dec_val=None):
    if dec_val is not None:
        cos_dec_err    *= -np.sin(dec_val)
        cov_ra_cos_dec *= -np.sin(dec_val)
    if (ra_err*cos_dec_err)**2 > cov_ra_cos_dec**2:
        return np.sqrt( (ra_err * cos_dec_err)**2 - cov_ra_cos_dec**2 ) * 2*PI * (180./PI)**2 * np.log(10)
    else:
        return None
