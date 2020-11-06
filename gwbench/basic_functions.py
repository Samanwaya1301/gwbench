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


from numpy import logical_and

#-----combine two sympy symbols strings without duplicates-----
def reduce_symbols_strings(string1,string2):
    # combine the two lists
    symbols_list = string1.split(' ') + string2.split(' ')
    # remove duplicates
    symbols_list = list(dict.fromkeys(symbols_list))

    # recreate a string of symbols and return it
    return ' '.join(symbols_list)

#-----delete symbols in sympy symbols strings, if not present in the other-----
def remove_symbols(string1,string2,keep_same=1):
    symbs_list1 = string1.split(' ')
    symbs_list2 = string2.split(' ')
    # remove unwanted symbols from 1
    if keep_same:
        symbols_list = [x for x in symbs_list1 if x in symbs_list2]
    else:
        symbols_list = [x for x in symbs_list1 if x not in symbs_list2]

    # recreate a string of symbols and return it
    return ' '.join(symbols_list)

#-----get subarray-----
def get_sub_array_ids(arr,sub_arr):
    return logical_and(arr>=sub_arr[0],arr<=sub_arr[-1])

#-----use only subset of dictionary-----
def get_sub_dict(dic,key_list,keep_in_list=1):
    if type(key_list) == str:
        key_list = key_list.split(' ')

    if keep_in_list:
        return {k:v for k,v in dic.items() if k in key_list}
    else:
        return {k:v for k,v in dic.items() if k not in key_list}
