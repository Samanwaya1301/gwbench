import numpy as np
from scipy.integrate import simps

#-----overlap function-------
def scalar_product_integrand(hf, gf, psd):
    temp = np.multiply(hf,np.conj(gf))
    return 2*np.divide(temp+np.conj(temp),psd)

def scalar_product_freq_array(hf, gf, psd, freqs, df=None):
    if df is None:
        return np.real(simps(scalar_product_integrand(hf, gf, psd),freqs))
    else:
        summ = 2.*np.real(((hf*np.conjugate(gf)+np.conjugate(hf)*gf)/psd).sum())*df
        return summ

#-----SNR function-------
def snr_square_integrand(hf, psd):
    return 4.*np.divide(np.power(np.abs(hf),2),psd)

def snr_square_freq_array(hf, psd, freqs):
    return simps(snr_square_integrand(hf, psd),freqs)

def snr_freq_array(hf, psd, freqs):
    return np.sqrt(snr_square_freq_array(hf, psd, freqs))

def snr_snr_sq_freq_array(hf, psd, freqs, df=None):
    if df is None:
        snr_sq = simps(snr_square_integrand(hf, psd),freqs)
        return np.sqrt(snr_sq), snr_sq
    else:
        summ = 4.*(np.divide(np.power(np.abs(hf),2),psd)).sum()*df
        summ = 2.*np.real(((hf*np.conjugate(hf)+np.conjugate(hf)*hf)/psd).sum())*df
        return np.sqrt(summ), summ

#-----fft method from Anuradha-------
def rfft_normalized(time_series, dt, n=None):
    return np.fft.rfft(time_series,n)*dt

def fft_normalized(time_series, dt, n=None):
    return np.fft.fft(time_series,n)*dt
