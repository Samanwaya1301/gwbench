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


#-----constants in SI units-----
GNewton = 6.6743e-11
cLight = 2.99792458e8
Msun = 1.9884099021470415e+30
Mpc = 3.085677581491367e+22
REarth = 6378136.6
AU = 1.4959787066e11
year = 3.1536e7

#-----convert mass in solar masses to seconds-----
MTsun = Msun * GNewton/cLight**3.
time_fac = MTsun
#-----convert mass/distance in solar mass/Mpc to dimensionless-----
strain_fac = GNewton/cLight**2.*Msun/Mpc
