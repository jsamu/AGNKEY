import numpy as np
from ccf_rm import *

DATA_DIR = './data'
DATA_INFO_FILE = 'main_data.txt'
flux_scale = 10.**(-15.)
bootstrap_iter = 10
######new params needed in new ccf_function 07/31/17
lag_range = 20
lag_step = 0.25

all_data = load_data(DATA_INFO_FILE, DATA_DIR)
line_data_list = all_data['LineName'].data
line_data_list = np.array(line_data_list)

peak = []
peak_err_high = []
peak_err_low = []
centroid = []
centroid_err_high =[]
centroid_err_low =[]
mass_bh = []


for obj in line_data_list:
	"""
	load in all info for an individual AGN
	"""
	obj_mask = all_data['LineName'] == obj
	obj_row = all_data[obj_mask]

	OBJ_NAME = obj
	PHOT_NAME = obj_row['PhotName'].data[0]
	flux_scaling_err = obj_row['ScaleErr'].data[0]
	print 'Working on %s'%PHOT_NAME

	hb_limits = (obj_row['HbetaWindow'].data[0]).rsplit('-')
	hb_limits = np.array([float(hb_limits[0]), float(hb_limits[1])])

	contin_limits1 = (obj_row['ContinWindow1'].data[0]).rsplit('-')
	contin_limits2 = (obj_row['ContinWindow2'].data[0]).rsplit('-')
	contin_limits = np.array([float(contin_limits1[0]), float(contin_limits1[1]), float(contin_limits2[0]), float(contin_limits2[1])])

	"""
	plot_spec: look at the first epoch's data
	line_disp: plot and fit RMS spectrum, extract Hbeta line dispersion
	light_curve: plot and return light curve data (V band, Hbeta)
	lag_dist: bootstrap lags, return average peak and centroid lags with 68% cl
	mass_BH: return the BH mass and errors
	"""
	agn_class = Agn(DATA_DIR, OBJ_NAME)
	#agn_class.plot_spec(hb_limits, contin_limits, flux_scale)
	disp = agn_class.line_disp(hb_limits, contin_limits, DATA_DIR, OBJ_NAME, flux_scale)
	v_date, v_mag, v_err, hb_date, hb_flux, hb_err = agn_class.lightcurve(hb_limits, contin_limits, DATA_DIR, OBJ_NAME, PHOT_NAME, flux_scale, flux_scaling_err)
	peak_avg, peak_errs, centroid_avg, cent_errs = lag_dist(bootstrap_iter, lag_range, lag_step, v_date, v_mag, v_err, hb_date, hb_flux, hb_err)
	m_bh = mass_BH(centroid_avg, disp)
	m_bh = '%.2e'%(m_bh)

	peak.append(peak_avg)
	peak_err_high.append(peak_errs[0])
	peak_err_low.append(peak_errs[1])
	centroid.append(centroid_avg)
	centroid_err_high.append(cent_errs[0])
	centroid_err_low.append(cent_errs[1])
	mass_bh.append(m_bh)


#create and output lag data table
lt = lag_table(all_data['PhotName'].data, peak, peak_err_high, peak_err_low, centroid, centroid_err_high, centroid_err_low, mass_bh)
print lt
