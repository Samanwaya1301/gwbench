## Examples:

**Beware that the number, passed to `num_gw_benchmarking.py`, `sym_gw_benchmarking.py`, and `multi_network.py`, needs to be smaller than `num_injs` as specified in each script.**

- Basic script to calculate antenna patterns, location phase factors, and PSDs:  
`python compute_ant_pat_lpf_psd.py`  

- Basic GW benchmarking with numeric derivatives:  
`python quick_start.py`  

- More elaborate GW benchmarking with numeric derivatives:  
`python num_gw_benchmarking.py 0`  

- Needed before symbolic GW benchmarking:  
`python generate_lambdified_functions.py`  

- More elaborate GW benchmarking with symbolic derivatives (otherwise the same as the numeric example):  
`python sym_gw_benchmarking.py 0`  

- Multi network GW benchmarking, setup with symbolic derivatives (can be switched to numeric):  
`python multi_network.py 0`  


**Available detector locations:**  
- standard sites:  
'H', 'L', 'V', 'K', 'I'

- fiducial sites:  
'ET1', 'ET2', 'ET3', 'U', 'A', 'W', 'B', 'C', 'N', 'S'  

**Available detector technologies (PSDs):**  
- 2G/2G+:  
'aLIGO', 'A+', 'V+', 'K+'  

- Voyager:  
'Voyager-CBO', 'Voyager-PMO'  

- 3G:  
'ET', 'CEwb'  
'CE1-10-CBO', 'CE1-20-CBO', 'CE1-30-CBO', 'CE1-40-CBO'  
'CE1-10-PMO', 'CE1-20-PMO', 'CE1-30-PMO', 'CE1-40-PMO'  
'CE2-10-CBO', 'CE2-20-CBO', 'CE2-30-CBO', 'CE2-40-CBO'  
'CE2-10-PMO', 'CE2-20-PMO', 'CE2-30-PMO', 'CE2-40-PMO'  

- LISA:  
'LISA-17', 'LISA-Babak17', 'LISA-Robson18'  

**DETECTOR KEYS TO PUT IN CODE**
- LIGO, VIRGO : ['aLIGO_H', 'aLIGO_L', 'aLIGO_V']
- ET : ['ET_ET1', 'ET_ET2', 'ET_ET3']
- CE : ['CE1-40-CBO_C', 'CE1-40-CBO_N', 'CE1-40-CBO_S'] or ['CE2-40-CBO_C', 'CE2-40-CBO_N', 'CE2-40-CBO_S']
