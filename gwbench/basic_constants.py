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
