import numpy as np
from astropy.io import ascii
from astropy.table import Table
#import seaborn as sns
from matplotlib import pyplot as pl
import scipy
from astropy.modeling import models, fitting
import collections

#sns.set_style("whitegrid")
#sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})


class Agn(object):
    """
    Load info from each epoch. Info file includes dates in HJD, seeing, and airmass corresponding to each spectrum.
    """
    def __init__(self, DATA_DIR, OBJ_NAME): 
        info = ascii.read('%s/lamp2008reduced/%s.scale/%s.scale.info.dat'%(DATA_DIR, OBJ_NAME, OBJ_NAME), delimiter='\t')
        self.filename = info['Filename'].data
        self.hjd = info['HJD at midpt'].data
        self.num_epochs = len(self.filename)
            
    def load_epochs(self, DATA_DIR, OBJ_NAME):
        """
        Load in all epochs of data as a dictionary with names of data files as the keys.
        Data files contain 3 columns: wavelength, flux density (in f-lambda units of 10^-15 erg/s/cm^2/A),
        and uncertainty in flux density. 
        """
        epoch_dict = {}
        
        for spec in self.filename:

            if spec not in epoch_dict:
                epoch_dict[spec] = ascii.read('%s/lamp2008reduced/%s.scale/%s'%(DATA_DIR,OBJ_NAME,spec))
            else:
                pass

        return epoch_dict
    
    def line_disp(self, hb_limits, bkgd_limits, DATA_DIR, OBJ_NAME, flux_scale):
        """
        Plots the spectrum of the first epoch as flux density vs. wavelength. Then fits
        a Gaussian to the Hbeta emission line and returns the velocity dispersion of the line.
        """
        epoch_dict_rmsd = self.load_epochs(DATA_DIR, OBJ_NAME)
        fluxdata_rmsd, flux_series_rmsd = self.stack_data(epoch_dict_rmsd, self.filename, flux_scale)
        wavelength = fluxdata_rmsd[0]

        #subtract linear continuum from spectrum of each epoch
        continuum_subtracted_fluxseries = []

        for epoch in flux_series_rmsd:
            continuum_fit = continuum_linfit(epoch, bkgd_limits)
            continuum_flux = continuum_fit[1] + continuum_fit[0]*epoch[0]
            continuum_subtracted_flux = epoch[1] - continuum_flux

            continuum_subtracted_epoch = np.vstack([epoch[0], continuum_subtracted_flux])
            continuum_subtracted_fluxseries.append(continuum_subtracted_epoch)

        continuum_subtracted_fluxseries = np.array(continuum_subtracted_fluxseries)

        #calculate mean and RMS deviation spectra
        mean_flux = np.average(continuum_subtracted_fluxseries[:,1,:], axis=0)
        mean_spectrum = np.vstack([wavelength, mean_flux])
        rmsd = np.sqrt((1./self.num_epochs)*np.sum((continuum_subtracted_fluxseries[:,1,:] - mean_flux)**2., axis=0))
        rmsd_spectrum = np.vstack([wavelength, rmsd])

        #constrain the spectra to only the line region with relevant surrounding continuum
        relevant_mask = (wavelength >= bkgd_limits[0]) & (wavelength <= bkgd_limits[3])
        wavelength = wavelength[relevant_mask]
        mean_spectrum = np.compress(relevant_mask, mean_spectrum, axis=1)
        rmsd_spectrum = np.compress(relevant_mask, rmsd_spectrum, axis=1)

        #select emission line region of data to be fit
        line_mask = (wavelength >= hb_limits[0]) & (wavelength <= hb_limits[1])
        wavelength_fit = wavelength[line_mask]
        flux_fit = rmsd_spectrum[1][line_mask]

        #fit gaussian to hbeta line
        #these parameters need to not be hard coded here
        mean_guess = np.median(wavelength_fit)
        amplitude_guess = np.max(flux_fit)
        stddev_guess = (hb_limits[1] - mean_guess)/2.

        g_init = models.Gaussian1D(amplitude=(amplitude_guess), mean=mean_guess, stddev=stddev_guess)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, wavelength_fit, flux_fit)

        fwhm = 2.*np.sqrt(2.*np.log(2.))*g.stddev.value
        c = 3.*10.**5.
        lambda0 = g.mean.value
        fwhm2 = fwhm*c/lambda0

        v_sig = g.stddev.value
        v_sig_km = v_sig*c/lambda0
        #I think we only care about line of sight velocity dispersion here
        #v_disp = v_sig_km*np.sqrt(3.)
        
        f, (ax1, ax2) = pl.subplots(2, 1, sharex=True, figsize=(11,7))

        ax1.plot(mean_spectrum[0], mean_spectrum[1])

        ax2.plot(rmsd_spectrum[0], rmsd_spectrum[1], 'b')
        ax2.plot(wavelength_fit, flux_fit, 'k')
        ax2.plot(wavelength_fit, g(wavelength_fit), 'r', label='Gaussian')
        #ax2.title('RMS spectrum fit with Gaussian (Object: %s)'%(OBJ_NAME))
        #ax2.xlabel('Wavelength ($\AA$)')
        #ax2.ylabel('Flux ($erg/s/cm^2$)')
        #ax2.text(min(rmsd_spectrum[0]), 0.8*max(flux_fit), r'$FWHM = %.2f \AA = %i  km/s$'%(fwhm, fwhm2))
        #ax2.text(min(rmsd_spectrum[0]), 0.6*max(flux_fit), r'$\sigma = %.2f \AA = %i  km/s$'%(v_sig, v_sig_km))

        #pl.savefig('./figures/lightcurves.png')
        pl.show()
        
        return v_sig_km

    def plot_spec(self, hb_limits, bkgd_limits, flux_scale):
        """
        Plots the spectrum of only the first epoch as flux density vs. wavelength.
        """
        epoch_dict1 = self.load_epochs()
        fluxdata1 = stack_data(epoch_dict1, self.filename[0], flux_scale)

        pl.figure(figsize=(12,8))
        pl.plot(fluxdata1[0], fluxdata1[1], 'k')
        pl.title(r'Sample spectrum from H$\beta$ time series')
        pl.xlabel('Wavelength ($\AA$)')
        pl.ylabel('Flux ($erg/s/cm^2$)')
        pl.savefig('./figures/singlespec.png')
        pl.show()
        
    def lightcurve(self, hb_limits, bkgd_limits, DATA_DIR, OBJ_NAME, PHOT_NAME, flux_scale, flux_scaling_err):
        """
        Returns the Hbeta continuum light curves as variables and plots them. Hbeta flux for each
        epoch is found by summing the spectrum within the hbmin to hbmax region, and subtracting a linear
        continuum fit within the window bmin to bmax. Continuum light curve is formed in band_lightcurve
        by taking data from 
        """
        v_date, v_mag, v_mag_err = band_lightcurve(DATA_DIR, OBJ_NAME, PHOT_NAME)

        epoch_dict = self.load_epochs(DATA_DIR, OBJ_NAME)
    
        hbsums = []
        lowlims = []
        uplims = []

        for epoch in self.filename:
            fluxdata = stack_data(epoch_dict, epoch, flux_scale)
            hbeta_fl = line_emission(fluxdata, hb_limits, bkgd_limits)
            lower_lim, upper_lim = flux_error(fluxdata, hb_limits, bkgd_limits)
            
            hbsums.append(hbeta_fl)
            lowlims.append(abs(hbeta_fl-lower_lim))
            uplims.append(abs(upper_lim-hbeta_fl))

        hb_sums = np.array(hbsums)
        lowlims = np.array(lowlims)
        uplims = np.array(uplims)
        limits = np.vstack([lowlims, uplims])
        hb_dates = np.array(self.hjd)
        hb_error = scale_error(hb_sums, lowlims, uplims, flux_scaling_err)
        
        
        
        f, (ax1, ax2) = pl.subplots(2, 1, sharex=True, figsize=(11,7))
        pl.subplots_adjust(hspace=0.07)

        ax1.errorbar(v_date, v_mag, yerr=v_mag_err, color='r', fmt='.')
        #ax1.invert_yaxis()
        ax1.set_ylabel(r'V band (mag)')

        ax2.errorbar(hb_dates,hb_sums, yerr=hb_error, color='b', fmt='.')
        ax2.set_ylabel(r'H$\beta$ Flux ($erg/s/cm^2$)')
        ax2.set_xlabel('HJD')

        #pl.savefig('./figures/lightcurves.png')
        pl.show()
        
        return v_date, v_mag, v_mag_err, hb_dates, hb_sums, hb_error

    def stack_data(self, epoch_dict, epoch, flux_scale):
        """
        Convert flux data table entries into a 3-dim array (wavelength, flux, fluxerror).
        This version is for the full data set.
        """
        data_table0 = epoch_dict[epoch[0]]
        #fluxdata = np.vstack((data_table0['col1'], data_table0['col1']*data_table0['col2']*flux_scale, data_table0['col1']*data_table0['col3']*flux_scale))
        #changed to keep data as flux density, instead of converting to flux
        fluxdata = np.vstack((data_table0['col1'], data_table0['col2']*flux_scale, data_table0['col3']*flux_scale))

        flux_series = np.empty([self.num_epochs, 3, fluxdata[0].size])
        for k, spec in enumerate(epoch):
            data_table = epoch_dict[spec]
            #fluxdata_ = np.vstack((data_table['col1'], data_table['col1']*data_table['col2']*flux_scale, data_table['col1']*data_table['col3']*flux_scale))
            fluxdata_ = np.vstack((data_table['col1'], data_table['col2']*flux_scale, data_table['col3']*flux_scale))
            flux_series[k] = fluxdata_

        return fluxdata, flux_series

def load_data(data_txt_file, DATA_DIR):
    """
    Loads in a text file of all AGN objects for use in the light curve function.
    """
    all_data_table = ascii.read('%s/%s'%(DATA_DIR, data_txt_file))

    return all_data_table

def stack_data(epoch_dict, epoch, flux_scale):
    """
    Convert flux data table entries into a 3-dim array (wavelength, flux, fluxerror).
    This version is for a single spectrum, used in plot_spec. Could be redundant.
    """
    data_table = epoch_dict[epoch]
    #fluxdata = np.vstack((data_table['col1'], data_table['col1']*data_table['col2']*flux_scale, data_table['col1']*data_table['col3']*flux_scale))
    #changed to keep data as flux density, instead of converting to flux
    fluxdata = np.vstack((data_table['col1'], data_table['col2']*flux_scale, data_table['col3']*flux_scale))

    return fluxdata

def line_emission(data, line_limits, bkgd_limits):
    """
    Finds emission line flux by summing total spectrum in emission line window, and subtracting
    sum of a linear continuum fit over same wavelength range. Takes as input 3 dimensional data
    array (wavelength, flux, fluxerror), wavelength window of the emission line as an array (min,max),
    and two wavelength windows for surrounding continuum regions as an array (min1,max1,min2,max2).
    """
    line = line_sum(data, line_limits)
    contin_fit = continuum_linfit(data, bkgd_limits)
    contin_sum = continuum_sum(data, line_limits, contin_fit)
    line_emission = line - contin_sum
    
    return line_emission

def line_sum(data, line_limits):
    """
    Returns a total sum underneath an emission line. Continuum included and needs
    to be subtracted later.
    """
    x = data[0]
    line_mask = (x >= line_limits[0]) & (x <= line_limits[1])
    line_comp = np.compress(line_mask,data,axis=1)
    line_total = np.trapz(line_comp[1],line_comp[0])
    return line_total

def continuum_linfit(data, bkgd_limits):
    """
    Returns parameters for a linear fit to a region given by bkgd_limits. The region is specified by
    4 wavelength values determining 2 wavelength windows: one to the left of the emission line and one
    to the right, trying to avoid any other lines/features.
    """
    x = data[0]
    contin_mask1 = (x >= bkgd_limits[0])&(x <= bkgd_limits[1])
    contin_mask2 =  (x >= bkgd_limits[2])&(x <= bkgd_limits[3])
    contin_comp = np.compress(contin_mask1+contin_mask2, data, axis=1)
    contin_fit = np.polyfit(contin_comp[0], contin_comp[1], 1)
    
    return contin_fit

def continuum_sum(data, line_limits, contin_fit):
    """
    Returns the sum under the continuum linear fit, but only in the region of the emission line_limits.
    """
    x = data[0]
    line_mask = (x >= line_limits[0]) & (x <= line_limits[1])
    line_comp = np.compress(line_mask,data,axis=1)
    contin_model = contin_fit[1] + contin_fit[0]*line_comp[0]
    contin_sum = np.trapz(contin_model,line_comp[0])

    return contin_sum

def flux_error(data, hb_limits, bkgd_limits):
    """
    Light curve error bars computed by forming a light curve for the lower and upper flux errors as given
    in LAMP data files. Returns an array of lower limits and an array of upper limits.
    """
    error_limits = []
    lower = np.vstack((data[0],data[1]-data[2],data[2]))
    upper = np.vstack((data[0],data[1]+data[2],data[2]))
    
    for i in [upper,lower]:
        data = i
        hbeta_fl = line_emission(data, hb_limits, bkgd_limits)
        error_limits.append(hbeta_fl)
        
    
    lower_lim = error_limits[0]
    upper_lim = error_limits[1]
    
    return lower_lim, upper_lim

def scale_error(hbsums, lowlims, uplims, sigma_nx):
    """
    An additional source of error in flux data from LAMP, which comes from the flux scaling
    using the OIII lines. See LAMP paper arXiv:1503.01146 for a detailed discussion of how
    sigma_nx is calculated. Previous lower and upper limits of flux are averaged here so only
    one final error array is returned. Could be a source of error itself, but likely not big.
    """
    avg_flux = np.mean(hbsums)/5.
    lims = np.vstack((lowlims, uplims))
    avg_lims = np.mean(lims, axis=0)
    scale_err = np.full(len(hbsums), sigma_nx*avg_flux)
    final_error = np.sqrt(avg_lims**2. + scale_err**2.)
    
    return final_error

def band_lightcurve(DATA_DIR, OBJ_NAME, PHOT_NAME):
    """
    Takes single LAMP 2008 photometry file in, and uses V band data. Any missing values are thrown out.
    Nights with more than one measurement are averaged and their errors weighted by n^(-1/2), where n
    is the number of epochs per night. Returns V band lightcurve with single point per night.
    Constant 2454000. days added to HJD from data file, see data file header.
    """
    photometry = ascii.read('%s/photometry.txt'%(DATA_DIR), format='cds', fill_values=[('--', np.NaN)])
    obj_mask = photometry['Name'] == PHOT_NAME
    photometry = photometry[obj_mask]
    photometry = photometry[~np.isnan(photometry['HJD-V'])]

    #maybe add an if statement to reject initial monitoring points
    if PHOT_NAME == 'Arp 151':
        v_date = np.array(photometry['HJD-V'].data) + 2454000.0
    else:
        v_date = np.array(photometry['HJD-V'].data) + 4000.0
    v_mag = np.array(photometry['Vmag'].data)
    v_mag_err = np.array(photometry['e_Vmag'].data)

    v_mag_avg_dict = {}
    for date in v_date:
        int_date = int(date)
        if int_date not in v_mag_avg_dict:
            v_mag_avg_dict[int_date] = [date]
        else:
            v_mag_avg_dict[int_date].append(date)

    #need this to counteract the (dis)order in which data was added to the dictionary
    v_mag_avg_dict = collections.OrderedDict(sorted(v_mag_avg_dict.items()))

    v_date_avg = []
    v_mag_avg = []
    v_mag_err_avg = []
    for key in v_mag_avg_dict.keys():
        night_dates = v_mag_avg_dict[key]
        if len(night_dates) == 1:
            indx = np.where(v_date == float(night_dates[0]))[0]
            v_date_avg.append(v_date[indx[0]])
            v_mag_avg.append(v_mag[indx[0]])
            v_mag_err_avg.append(v_mag_err[indx[0]])

        elif len(night_dates) > 1:
            v_date_avg_p = []
            v_mag_avg_p = []
            v_mag_err_avg_p = []

            for i in night_dates:
                indx = np.where(v_date == float(i))[0]
                v_date_avg_p.append(v_date[indx[0]])
                v_mag_avg_p.append(v_mag[indx[0]])
                v_mag_err_avg_p.append(v_mag_err[indx[0]])

            #print v_date_avg_p, v_mag_avg_p, v_mag_err_avg_p
            v_date_avg_p = np.array(v_date_avg_p)
            v_mag_avg_p = np.array(v_mag_avg_p)
            v_mag_err_avg_p = np.array(v_mag_err_avg_p)

            v_date_avg.append(np.average(v_date_avg_p))
            v_mag_avg.append(np.average(v_mag_avg_p))
            v_mag_err_avg.append(np.average(v_mag_err_avg_p)/np.sqrt(v_date_avg_p.size))

    v_date_avg = np.array(v_date_avg)
    v_mag_avg = np.array(v_mag_avg)
    v_mag_err_avg = np.array(v_mag_err_avg)

    #convert from magnitude to flux density
    #V_band_zero_point = 3.631*10.**(-9.)#erg/s/cm^2/A from Bessel 1998
    V_band_zero_point = 3.735*10.**(-9.)#erg/s/cm^2/A from Mann 2015
    zp_lambda_eff = 5529.#angstroms
    V_band_flux = zp_lambda_eff*V_band_zero_point*10.**(-v_mag_avg/2.5)
    V_band_flux_with_err = zp_lambda_eff*V_band_zero_point*10.**(-(v_mag_avg-v_mag_err_avg)/2.5)
    V_band_flux_err = V_band_flux_with_err - V_band_flux

    #return v_date_avg, v_mag_avg, v_mag_err_avg
    return v_date_avg, V_band_flux, V_band_flux_err

def ccfprep_IC(x_a, a, a_err, x_b, b, b_err):
    """
    Interpolates continuum light curve ("a") to match the average sampling of the shifted line
    light curve ("b"). Returns both light curves with only the continuum being interpolated.
    """
    #negative accounts for the effect of "a" being a magnitude
    #a = -a
    x_a = x_a - 2454520.
    x_b = x_b - 2454520.

    #match the exact sampling of the line light curve in region of overlap
    time_overlap_mask1 = x_b >= x_a[0]
    time_overlap_mask2 = x_b <= x_a[-1]
    time_overlap_mask = time_overlap_mask1 & time_overlap_mask2
    x_a_interp = x_b[time_overlap_mask]

    a_interp = np.interp(x_a_interp,x_a,a)
    a_err_interp = np.interp(x_a_interp, x_a, a_err)
    
    return x_a_interp, a_interp, a_err_interp, x_b[time_overlap_mask], b[time_overlap_mask], b_err[time_overlap_mask]

def ccfprep_IL(x_a, a, a_err, x_b, b, b_err):
    """
    Interpolates shifted line emisson light curve ("b") to match the average sampling of the
    continuum light curve ("a"). Returns both light curves with only the line one
    being interpolated.
    """
    #negative accounts for the effect of "a" being a magnitude
    #a = -a
    x_a = x_a - 2454520.
    x_b = x_b - 2454520.

    #match the exact sampling of the continuum light curve
    time_overlap_mask1 = x_a >= x_b[0]
    time_overlap_mask2 = x_a <= x_b[-1]
    time_overlap_mask = time_overlap_mask1 & time_overlap_mask2
    x_b_interp = x_a[time_overlap_mask]

    b_interp = np.interp(x_b_interp,x_b,b)
    b_err_interp = np.interp(x_b_interp, x_b, b_err)
    
    return x_a[time_overlap_mask], a[time_overlap_mask], a_err[time_overlap_mask], x_b_interp, b_interp, b_err_interp

def ccf_sum(y1, y1_err, y1_std, y2, y2_err, y2_std):
    """
    Calculates the CCF for an individual lag. Inputs are determined in ccf_full by how much
    of each time series is overalpping for the particular lag. Returns a single scalar value
    which is appended to an array of CCF values in ccf_full.
    """
    N_ccf = y1.size
    y1_meansub = y1 - np.average(y1)
    y2_meansub = y2 - np.average(y2)
    ccf_w = (1./N_ccf)*np.sum( (y1_meansub*y2_meansub)/(y1_std*y2_std) )

    return ccf_w

def ccf_new(lag_range, lag_step, v_date_sample, v_sample, v_err_sample, hb_date_sample, hb_sample, hb_err_sample):

    ccf_IC = []
    ccf_IL = []

    lags = np.arange(-lag_range, lag_range + lag_step, lag_step)

    for i in lags:
        #shift line emission light curve by each lag
        hb_date_sample_shifted = hb_date_sample + i

        #interpolate continuum light curve to match the timing of the line emission light curve
        v_dateIC, v_magIC, v_errIC, hb_dateIC, hb_fluxIC, hb_errIC = ccfprep_IC(v_date_sample, v_sample, v_err_sample, hb_date_sample_shifted, hb_sample, hb_err_sample)
        
        #interpolate line emission light curve to match the timing of the continuum light curve
        v_dateIL, v_magIL, v_errIL, hb_dateIL, hb_fluxIL, hb_errIL = ccfprep_IL(v_date_sample, v_sample, v_err_sample, hb_date_sample_shifted, hb_sample, hb_err_sample)

        #make sure at least 2 data points are overlapping for CCF calculationg
        if (len(v_dateIC) <= 1) or (len(v_dateIL) <= 1):
            print 'Too few overlapping points to calculate CCF'
            ccf_IC.append(0.0)
            ccf_IL.append(0.0)
            continue
        else:
            ccf_IC_value = ccf_sum(v_magIC, v_errIC, np.std(v_magIC), hb_fluxIC, hb_errIC, np.std(hb_fluxIC))
            ccf_IC.append(ccf_IC_value)

            ccf_IL_value = ccf_sum(v_magIL, v_errIL, np.std(v_magIL), hb_fluxIL, hb_errIL, np.std(hb_fluxIL))
            ccf_IL.append(ccf_IL_value)

    ccf_IC = np.array(ccf_IC)
    ccf_IL = np.array(ccf_IL)

    ccf_avg_new = np.average([ccf_IC, ccf_IL], axis=0)

    '''
    plot needs to be edited, copied and pasted from old ccf function
    f, ((ax1, ax3)) = pl.subplots(2, 1, figsize=(7,12))
    f.subplots_adjust(hspace=0.3)

    ax1.set_title('Interpolated Time Series')
    ax1.set_xlabel(r'Time (HJD)')
    ax1.plot(x1,y1_mean_sub, 'k', label='$C(t)$')
    #ax1.plot(x1,y1_mean_sub, 'r--', label='$C(t)$')
    ax1.plot(x2,y2_mean_sub, 'b', label=r'$L(t+\tau_{true})$')
    ax1.plot(x2 + maxlag_norm, y2_mean_sub, 'r--', label=r'$L(t+\tau_{peak})$')
    ax1.legend(loc=2)
    ax1.set_ylabel('Difference from mean')

    ax3.set_ylabel('Normalized correlation coefficient')
    ax3.set_title(r'CCF($\tau$)')
    ax3.set_xlabel(r'$\tau$ (days)')
    ax3.plot(time_limit, ccf_norm_limit)
    ax3.text(0.25*max(time_limit), 0.9*max(ccf_norm_limit), r'$\tau_{peak}=%.3f$'%(maxlag_norm))
    
    ax4.set_title('Un-shifted time series', fontsize=20)
    ax4.invert_yaxis()
    ax4.tick_params(labelsize=16)
    ax4.plot(x1, y1_mean_sub, label=r'$C(t)$')
    ax4.plot(x2 + maxlag_norm, y2_mean_sub, '--', label=r'$C(t-\tau_{peak})$')
    ax4.legend(loc=2, fontsize=18)
    
    #pl.savefig('./figures/correlation.png')    
    pl.show()
    '''

    return lags, ccf_avg_new

def mass_BH(lag, dispersion):
    G_const = 6.674*10.**(-20.)#km^3/(kgs^2)
    c = 2.99*10.**5.
    lightdays2sec = 86400.
    m_BH = 5.5*(c*lag*lightdays2sec*dispersion**2.)/G_const
    m_sun = 1.988*10.**30.
    m_BH_solar = m_BH/m_sun
    #m_BH_solar = float("{0:.2e}".format(m_BH_solar))
    #print 'M_bh = %.2e M_sun'%(m_BH_solar)

    return m_BH_solar

def bootstrap_resample(x, y, y_err, n=None):
    """ 
    Resample the time series y with replacement. If a point is chosen more than
    once, its error is weighted by n^(-1/2) where n is the number of times it was chosen.
    Also adds in "flux randomization" where a random Gaussian deviation is added to the
    flux value (see arXiv:0908.0003).
    """
    if n == None:
        n = len(x)

    resample_i = sorted(np.random.choice(n, n))
    x_resample = list(x[resample_i])

    sample_counter = np.empty(n)
    for i,j in enumerate(x):
        sample_counter[i] = x_resample.count(j)

    sample_mask = sample_counter > 0
    masked_sample_counter = sample_counter[sample_mask]
    x_resample = x[sample_mask]
    y_err_resample = y_err[sample_mask]/np.sqrt(masked_sample_counter)
    y_resample = y[sample_mask] + np.random.normal(0.,1., y_err_resample.size)*y_err_resample

    return x_resample, y_resample, y_err_resample

def lag_dist(n_inter, lag_range, lag_step, v_date, v, v_err, hb_date, hb_flux, hb_err):
    """
    Calculate the CCF over n_iter iterations (for resampled light curves) twice: once while interpolating
    continuum and once while interpolating line emission. Returns means of the distributions
    for peak lag and centroid lag, as well as their 68% confidence intervals. Plots distributions.
    """
    num_iter = np.arange(n_inter)

    #new peak and centroid lags from revised CCF 07/24/17
    peak_lags = np.empty(n_inter)
    cent_lags = np.empty(n_inter)

    for i in num_iter:
        #resample each light curve first
        v_date_sample, v_sample, v_err_sample = bootstrap_resample(v_date, v, v_err)
        hb_date_sample, hb_sample, hb_err_sample = bootstrap_resample(hb_date, hb_flux, hb_err)

        #calculate ccf
        lags, ccf = ccf_new(lag_range, lag_step, v_date_sample, v_sample, v_err_sample, hb_date_sample, hb_sample, hb_err_sample)

        #######################07/24/17
        #find peak lag
        peak_lags[i] = lags[np.argmax(ccf)]

        #find  and centroid lag
        max_r_mask = ccf >= 0.8*np.max(ccf)
        #print max_r_mask
        ccf_max_r = ccf[max_r_mask]
        lags_max_r = lags[max_r_mask]
        centroid_numerator = np.trapz(lags_max_r*ccf_max_r, lags_max_r)
        centroid_denominator = np.trapz(ccf_max_r, lags_max_r)
        centroid = centroid_numerator/centroid_denominator
        centroid = float('%.3f'%centroid)
        cent_lags[i] = centroid

    #some centroids have strange behavior(?) so NaN's are masked out of all centroid arrays
    #cent_lags = cent_lags[~np.isnan(cent_lags)]

    #find means of new distributions 07/24/17
    peak_lags_mean = np.average(peak_lags)
    cent_lags_mean = np.average(cent_lags)

    #calculate new confidence intervals 07/24/17
    peak_lags.sort()
    p_e1 = peak_lags_mean - peak_lags[int(0.16*(peak_lags.size-1))]
    p_e2 = peak_lags_mean - peak_lags[int(0.84*(peak_lags.size-1))]
    peak_errors = [float("{0:.2f}".format(-p_e2)), float("{0:.2f}".format(p_e1))]

    cent_lags.sort()
    c_e1 = cent_lags_mean - cent_lags[int(0.16*(cent_lags.size-1))]
    c_e2 = cent_lags_mean - cent_lags[int(0.84*(cent_lags.size-1))]
    cent_errors = [float("{0:.2f}".format(-c_e2)), float("{0:.2f}".format(c_e1))]


    '''
    #plot distributions of lags found after averaging the ccf's
    fig2, (ax21, ax22) = pl.subplots(2, 1, figsize=(7,10))
    pl.subplots_adjust(hspace=0.4)
    ax21.set_title(r'$\bar{\tau}_{peak}=%.2f\binom{+%.2f}{-%.2f}$'%(peak_lags_mean, abs(p_e2), abs(p_e1)))
    ax21.set_xlabel('Lag (days)')
    ax21.set_ylabel('Peak frequency')
    sns.distplot(peak_lags, kde=False, bins=30, color='b', ax=ax21, hist_kws={"linewidth": 1.0, "alpha": 0.75})

    ax22.set_xlabel('Lag (days)')
    ax22.set_ylabel('Centroid frequency')
    ax22.set_title(r'$\bar{\tau}_{cent}=%.2f\binom{+%.2f}{-%.2f}$'%(cent_lags_mean, abs(c_e2), abs(c_e1)))
    sns.distplot(cent_lags, kde=False, bins=30, color='b', ax=ax22, hist_kws={"linewidth": 1.0, "alpha": 0.75})

    #pl.savefig('./figures/bootstrap_avccf.png')
    pl.show()
    '''

    return peak_lags_mean, peak_errors, cent_lags_mean, cent_errors

def lag_table(PHOT_NAME, peak, peak_err_high, peak_err_low, centroid, centroid_err_high, centroid_err_low, mass_bh):
    lag_data_table = Table([PHOT_NAME, peak, peak_err_high, peak_err_low, centroid, centroid_err_high, centroid_err_low, mass_bh],
                         names=['Name', 'PeakLag', 'PErr(up)', 'PErr(low)', 'CentLag', 'CErr(up)', 'CErr(low)', 'M_BH'])
    ascii.write(lag_data_table, 'lag_data.txt', delimiter='\t')

    return lag_data_table
