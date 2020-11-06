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


from copy import copy

import sympy as sp

from gwbench.wf_models import lal_bbh_np
from gwbench.wf_models import lal_bns_np
from gwbench.wf_models import tf2_np
from gwbench.wf_models import tf2_sp
from gwbench.wf_models import tf2_tidal_np
from gwbench.wf_models import tf2_tidal_sp

class Waveform(object):

    ###
    #-----Init methods-----
    def __init__(self, wf_model_name=None, wf_other_var_dic=None):
        if wf_model_name == None:
            wf_symbs_string = None
            hfpc_np = None
            hfpc_sp = None
        else:
            wf_symbs_string, hfpc_np, hfpc_sp = select_wf_model_quants(wf_model_name)

        self.wf_model_name = wf_model_name
        self.wf_other_var_dic = wf_other_var_dic
        self.wf_symbs_string = wf_symbs_string
        self.hfpc_np = hfpc_np
        self.hfpc_sp = hfpc_sp


    ###
    #-----Getter methods-----
    def get_sp_expr(self):
        symb_dic = {}
        for name in self.wf_symbs_string.split(' '):
            symb_dic[name] = sp.symbols(name,real=True)

        if self.wf_other_var_dic == None:
            return self.hfpc_sp(*list(symb_dic.values()))
        else:
            return self.hfpc_sp(*list(symb_dic.values()),*list(self.wf_other_var_dic.values()))

    def eval_np_func(self,f,inj_params):
        if isinstance(inj_params, dict):
            if self.wf_other_var_dic is None:
                return self.hfpc_np(f,**inj_params)
            else:
                return self.hfpc_np(f,**inj_params,**self.wf_other_var_dic)
        elif isinstance(inj_params, list):
            if self.wf_other_var_dic is None:
                return self.hfpc_np(f,*inj_params)
            else:
                return self.hfpc_np(f,*inj_params,**self.wf_other_var_dic)


    ###
    #-----IO methods-----
    def print_waveform(self):
        for key,value in vars(self).items():
            print(key.ljust(16,' '),'  ',value)
            print()


###
#-----Get waveform functions for np, sp and the symbols string based on the model name-----
def select_wf_model_quants(wf_model_name):
    if wf_model_name == 'lal_bbh':
        np_mod = lal_bbh_np
        sp_mod = None
    elif wf_model_name == 'lal_bns':
        np_mod = lal_bns_np
        sp_mod = None
    elif wf_model_name == 'tf2':
        np_mod = tf2_np
        sp_mod = tf2_sp
    elif wf_model_name == 'tf2_tidal':
        np_mod = tf2_tidal_np
        sp_mod = tf2_tidal_sp

    if sp_mod is None:
        return np_mod.wf_symbs_string, np_mod.hfpc, None
    elif np_mod is None:
        return sp_mod.wf_symbs_string, None, sp_mod.hfpc
    else:
        return np_mod.wf_symbs_string, np_mod.hfpc, sp_mod.hfpc
