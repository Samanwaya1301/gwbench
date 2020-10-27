import os
from copy import copy

import astropy.cosmology as apcosm
import numpy as np
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

import gwbench.basic_relations as brs

PI = np.pi

###
#-----standard values for reference-----
std_vals = [{'source':'BBH', 'dist':'power_uniform', 'mmin':5, 'mmax':100, 'alpha':-1.6},
            {'source':'BNS', 'dist':'gaussian', 'mean':1.35, 'sigma':0.15},
            {'source':'BNS', 'dist':'uniform', 'mmin':0.8, 'mmax':3}]

###
#-----CBC parameters sampler-----
def injections_CBC_params_redshift(cosmo_dict,mass_dict,spin_dict,redshifted,num_injs=10,seed=None,file_path=None):
    rng = np.random.default_rng(seed)
    seeds = rng.integers(100000,size=4)

    m1_vec, m2_vec = mass_sampler(mass_dict,num_injs,seeds[0])
    chi1x_vec, chi1y_vec, chi1z_vec, chi2x_vec, chi2y_vec, chi2z_vec = spin_sampler(spin_dict,num_injs,seeds[1])
    z_vec, DL_vec = redshift_lum_distance_sampler(cosmo_dict,num_injs,seeds[2])
    iota_vec, ra_vec, dec_vec, psi_vec = angle_sampler(num_injs,seeds[3])

    Mc_vec, eta_vec = get_Mc_eta(m1_vec,m2_vec)
    if redshifted: Mc_vec *= (1. + z_vec)

    params = [Mc_vec, eta_vec, chi1x_vec, chi1y_vec, chi1z_vec, chi2x_vec, chi2y_vec, chi2z_vec, DL_vec, iota_vec, ra_vec, dec_vec, psi_vec, z_vec]
    if file_path is not None: save_injections(params,file_path)
    return params


###
#-----IO functions-----
def load_injections(file_path):
    return np.transpose(np.loadtxt(file_path))

def save_injections(params,file_path):
    np.savetxt(os.path.join(file_path), np.transpose(np.array([params][0])), delimiter = ' ')


###
#-----angle samplers-----
def angle_sampler(num_injs,seed):
    rng = np.random.default_rng(seed)
    iota_vec = np.arccos(rng.uniform(low=-1, high=1, size=num_injs))
    ra_vec   = rng.uniform(low=0., high=2.*PI, size=num_injs)
    dec_vec  = np.arccos(rng.uniform(low=-1, high=1, size=num_injs)) - PI/2.
    psi_vec  = rng.uniform(low=0., high=2.*PI, size=num_injs)
    return iota_vec, ra_vec, dec_vec, psi_vec

###
#-----spin samplers-----
def spin_sampler(spin_dict,num_injs,seed):
    rng = np.random.default_rng(seed)
    chi_lo   = spin_dict['chi_lo']
    chi_hi   = spin_dict['chi_hi']
    dim = spin_dict['dim']

    if   dim == 1:
        chiz_vecs = [rng.uniform(low=chi_lo, high=chi_hi, size=num_injs) for i in range(2)]
        return [np.zeros(num_injs), np.zeros(num_injs), chiz_vecs[0], np.zeros(num_injs), np.zeros(num_injs), chiz_vecs[1]]
    elif dim == 3:
        if   spin_dict['geom'] == 'cartesian':
            return [rng.uniform(low=chi_lo, high=chi_hi, size=num_injs) for i in range(6)]
        elif spin_dict['geom'] == 'spherical':
            # chi1
            chi_vec   = (rng.uniform(low=chi_lo**3., high=chi_hi**3., size=num_injs))**(1./3.)
            theta_vec = np.arccos(rng.uniform(low=-1, high=1, size=num_injs))
            phi_vec   = rng.uniform(low=0., high=2.*PI, size=num_injs)
            chi1x_vec, chi1y_vec, chi1z_vec = get_cartesian_from_spherical(chi_vec,theta_vec,phi_vec)
            # chi2
            chi_vec   = (rng.uniform(low=chi_lo**3., high=chi_hi**3., size=num_injs))**(1./3.)
            theta_vec = np.arccos(rng.uniform(low=-1, high=1, size=num_injs))
            phi_vec   = rng.uniform(low=0., high=2.*PI, size=num_injs)
            chi2x_vec, chi2y_vec, chi2z_vec = get_cartesian_from_spherical(chi_vec,theta_vec,phi_vec)
            return [chi1x_vec, chi1y_vec, chi1z_vec, chi2x_vec, chi2y_vec, chi2z_vec]

def get_cartesian_from_spherical(r,theta,phi):
    return r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)


###
#-----mass samplers-----
def mass_sampler(mass_dict,num_injs,seed):
    rng = np.random.default_rng(seed)
    if mass_dict['dist'] == 'gaussian':
        mmin  = mass_dict['mmin']
        mmax  = mass_dict['mmax']
        mean  = mass_dict['mean']
        sigma = mass_dict['sigma']
        m1_m2 = 0

        m1_vec = np.zeros(num_injs)
        ids = np.arange(num_injs)
        while ids.size > 0:
            m1_vec[ids] = rng.normal(loc=mean, scale=sigma, size=ids.size)
            ids = np.nonzero(np.logical_not(np.logical_and(m1_vec > mmin, m1_vec < mmax)))[0]

        m2_vec = np.zeros(num_injs)
        ids = np.arange(num_injs)
        while ids.size > 0:
            m2_vec[ids] = rng.normal(loc=mean, scale=sigma, size=ids.size)
            ids = np.nonzero(np.logical_not(np.logical_and(m2_vec > mmin, m2_vec < mmax)))[0]

    elif mass_dict['dist'] == 'power':
        mmin  = mass_dict['mmin']
        mmax  = mass_dict['mmax']
        alpha = mass_dict['alpha'] + 1
        m1_m2 = 0

        m1_vec = (mmin**alpha + (mmax**alpha - mmin**alpha)*rng.random(num_injs))**(1./alpha)
        m2_vec = (mmin**alpha + (mmax**alpha - mmin**alpha)*rng.random(num_injs))**(1./alpha)

    elif mass_dict['dist'] == 'power_uniform':
        mmin  = mass_dict['mmin']
        mmax  = mass_dict['mmax']
        alpha = mass_dict['alpha'] + 1
        m1_m2 = 1

        m1_vec = (mmin**alpha + (mmax**alpha - mmin**alpha)*rng.random(num_injs))**(1./alpha)
        m2_vec = rng.uniform(mmin,m1_vec)

    elif mass_dict['dist'] == 'uniform':
        mmin  = mass_dict['mmin']
        mmax  = mass_dict['mmax']
        m1_m2 = 0

        m1_vec = rng.uniform(low=mmin, high=mmax, size=num_injs)
        m2_vec = rng.uniform(low=mmin, high=mmax, size=num_injs)

    if m1_m2: return m1_vec, m2_vec
    else:     return make_m1_m2(m1_vec,m2_vec,num_injs!=1)

#-----mass handling functions-----
def get_Mc_eta(m1_vec,m2_vec):
    eta_vec  = brs.eta_of_q(m1_vec/m2_vec)
    Mc_vec   = brs.Mc_of_M_eta(m1_vec+m2_vec,eta_vec)
    return Mc_vec, eta_vec

def make_m1_m2(m1,m2,vec=1):
    if vec:
        mt = copy(m1)
        ids = np.where(m1<m2)
        m1[ids] = m2[ids]
        m2[ids] = mt[ids]
        return m1, m2
    else:
        if m1 < m2:
            return m2, m1
        else:
            return m1, m2


###
#-----redshift and lum distance samplers-----
def redshift_lum_distance_sampler(cosmo_dict,num_injs,seed):
    zmin    = cosmo_dict['zmin']
    zmax    = cosmo_dict['zmax']

    keys = list(cosmo_dict.keys())
    if 'Om0' in keys: Om0 = cosmo_dict['Om0']
    else:                 Om0 = None
    if 'Ode0' in keys: Ode0 = cosmo_dict['Ode0']
    else:                 Ode0 = None
    if 'H0' in keys:      H0 = cosmo_dict['H0']
    else:                 H0 = None

    if None in (Om0,Ode0,H0):
        cosmo = apcosm.Planck15
    else:
        cosmo = apcosm.LambdaCDM(H0=H0, Om0=Om0, Ode0=Ode0)

    if cosmo_dict['sampler'] == 'inversion':
        nzs = None
        z_vec = uniform_comoving_volume_redshift_inversion_sampler(zmin,zmax,cosmo,num_injs,seed,nzs)
    elif cosmo_dict['sampler'] == 'rejection':
        nzs = 40
        z_vec = uniform_comoving_volume_redshift_rejection_sampler(zmin,zmax,cosmo,num_injs,seed,nzs)

    return z_vec, cosmo.luminosity_distance(z_vec).value

#-----redshift samplers-----
def uniform_comoving_volume_redshift_inversion_sampler(zmin,zmax,cosmo,num_injs,seed=None,nzs=None):
    rng = np.random.default_rng(seed)
    if nzs is None: nzs = max(5001, 50 * int((zmax-zmin)) + 1)
    zs = np.linspace(zmin,zmax,nzs)
    dist = (lambda z: ((4.*PI*cosmo.differential_comoving_volume(z).value)/(1.+z)))(zs)
    window_cdf = np.array([simps(dist[:i],zs[:i]) for i in range(1,nzs+1)]) / simps(dist,zs)
    inv_window_cdf = interp1d(window_cdf, zs)
    return inv_window_cdf(rng.random(num_injs))

def uniform_comoving_volume_redshift_rejection_sampler(zmin,zmax,cosmo,num_injs,seed=None,nzs=None):
    rng = np.random.default_rng(seed)
    dist = lambda z: ((4.*PI*cosmo.differential_comoving_volume(z).value)/(1.+z))
    window_norm = quad(dist, zmin, zmax)[0]
    flip_window_pdf = lambda z: -dist(z) / window_norm
    window_pdf_max = -minimize_scalar(flip_window_pdf,bounds=[zmin,zmax],method='bounded').fun

    if nzs is None or nzs < 2: nzs = 2
    zs = np.linspace(zmin,zmax,nzs)
    segment_nums = np.asarray((num_injs * np.array([-quad(flip_window_pdf,zs[i],zs[i+1])[0] for i in range(nzs-1)])), dtype=int)
    segment_pts = np.sum(segment_nums)
    segment_maxs = np.array([-minimize_scalar(flip_window_pdf,bounds=[zs[i],zs[i+1]],method='bounded').fun for i in range(nzs-1)])

    z_sample = np.zeros(num_injs)

    for j,num in enumerate(segment_nums):
        id_shift = np.sum(segment_nums[:j+1]) - num
        ids = np.arange(num) + id_shift
        while ids.size > 0:
            z_sample[ids] = rng.uniform(zs[j],zs[j+1],ids.size)
            ids = ids[np.nonzero(rng.uniform(0.,segment_maxs[j],ids.size) >= -flip_window_pdf(z_sample[ids]))[0]]

    if nzs > 2:
        ids = np.arange(num_injs - segment_pts) + segment_pts
        while ids.size > 0:
            z_sample[ids] = rng.uniform(zmin,zmax,ids.size)
            ids = ids[np.nonzero(rng.uniform(0.,window_pdf_max,ids.size) >= -flip_window_pdf(z_sample[ids]))[0]]
        rng.shuffle(z_sample)

    return z_sample
