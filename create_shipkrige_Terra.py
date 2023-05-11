"""
Create input data for kriging from CERES SSF1deg Terra data

Code for making kriging input data for 'Detection of large-scale cloud microphysical changes and evidence for decreasing cloud brightness within a major shipping corridor due to the 2020 International Maritime Organization marine fuel sulfur regulations'

Modifications
-------------
11 April 2023: Michael Diamond, Tallahassee, FL
    -Created
10 May 2023: Michael Diamond, Tallahassee, FL
    -Final form for ACP Letters initial submission
"""

#Import libraries
import numpy as np
import xarray as xr
from scipy import stats
from glob import glob
import os
import warnings

#Set paths
dir_data = '/Users/michaeldiamond/Documents/Data/'


"""
Create shipkrige nc file for R kriging code
"""

sk = xr.Dataset()

#Load legacy shipkrige file from Diamond et al. (2020), AGU Adv.
orig = xr.open_dataset('/Users/michaeldiamond/Dropbox/ShipTrackAnalysis_MichaelDiamond/origData/shipkrige_combinedvars.nc')
sk['month'] = orig['month']
sk['lat'] = orig['lat']
sk['lon'] = orig['lon']
sk['EDGAR_SO2'] = orig['EDGAR_SO2']

#
###Load and manipulate SSF1deg data
#

#Load data
Terra = xr.open_dataset(glob(dir_data+'CERES/SSF1deg/CERES_SSF1deg-Month_Terra-MODIS_Ed4.1_Subset_*.nc')[0])

#Calculate albedo, total cloud fraction, and overcast albedo
Terra['A'] = Terra['toa_sw_all_mon']/Terra['toa_solar_all_mon']

Terra['C'] = (Terra['cldarea_high_day_mon']+Terra['cldarea_mid_high_day_mon']+Terra['cldarea_mid_low_day_mon']+Terra['cldarea_low_day_mon'])/100

Terra['Acld'] = (Terra['A']-(1-Terra['C'])*.1)/Terra['C'] #Assumes Aclr = 0.1

Terra['Clow'] = Terra['cldarea_low_day_mon']/100


#
###Fill in data
#
print('Starting...')
for sat, ds in zip(['Terra'],[Terra]):
    print('~%s~' % sat)
    with warnings.catch_warnings(): #Ignore runtime warnings
        warnings.simplefilter("ignore")
    
        #
        ###Get "climatological" (2002-2019) means and remap grid
        #
        print('...clim...')

        sk['%s_cer_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_cer_clim' % sat].attrs = {'long_name' : '%s liquid cloud effective radius (3.7 µm) averaged from 2003 to 2019' % sat,'units' : 'µm'}

        sk['%s_cot_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_cot_clim' % sat].attrs = {'long_name' : '%s cloud optical thickness averaged from 2003 to 2019' % sat,'units' : '1'}

        sk['%s_lwp_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_lwp_clim' % sat].attrs = {'long_name' : '%s liquid water path (3.7 µm) averaged from 2003 to 2019' % sat,'units' : 'g m-2'}

        sk['%s_Ctot_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_Ctot_clim' % sat].attrs = {'long_name' : '%s total cloud fraction averaged from 2003 to 2019' % sat,'units' : '1'}

        sk['%s_Clow_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_Clow_clim' % sat].attrs = {'long_name' : '%s low cloud fraction averaged from 2003 to 2019' % sat,'units' : '1'}

        sk['%s_A_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_A_clim' % sat].attrs = {'long_name' : '%s all-sky albedo averaged from 2003 to 2019' % sat,'units' : '1'}

        sk['%s_Acld_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_Acld_clim' % sat].attrs = {'long_name' : '%s overcast albedo (assuming a clear-sky albedo of 0.1) averaged from 2003 to 2019' % sat,'units' : '1'}

        sk['%s_Tskin_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_Tskin_clim' % sat].attrs = {'long_name' : '%s skin temperature averaged from 2003 to 2019' % sat,'units' : 'K'}

        sk['%s_EIS_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_EIS_clim' % sat].attrs = {'long_name' : '%s Estimated Inversion Strength averaged from 2003 to 2019' % sat,'units' : 'K'}

        sk['%s_WS_clim' % sat] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
        sk['%s_WS_clim' % sat].attrs = {'long_name' : '%s wind speed averaged from 2003 to 2019' % sat,'units' : 'm s-1'}


        #Loop through months for temporal averaging

        ymask = np.logical_and(ds.time.dt.year>=2002,ds.time.dt.year<=2019)

        for m in range(1,13):

            mmask = ds.time.dt.month==m
            tmask = np.logical_and(ymask,mmask)

            for var1, var2 in zip(['cer','cot','lwp','Ctot','Clow','A','Acld','Tskin','EIS','WS'],['cldwatrad37_total_day_mon','cldtau_total_day_mon','lwp37_total_day_mon','C','Clow','A','Acld','aux_skint_mon','aux_inversion_mon','aux_wind_speed_mon']):
                sk['%s_%s_clim' % (sat,var1)][m-1,:,:180] = np.nanmean(ds[var2][tmask,::-1,:][:,:,ds.lon>180],axis=0)
                sk['%s_%s_clim' % (sat,var1)][m-1,:,180:] = np.nanmean(ds[var2][tmask,::-1,:][:,:,ds.lon<180],axis=0)


        #
        ###Take 3-year monthly averages and remap grid
        #

        for y in np.arange(2002,2023,3):
            print('...%s...' % y)

            sk['%s_cer_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_cer_%s' % (sat,y)].attrs = {'long_name' : '%s liquid cloud effective radius (3.7 µm) averaged from %s to %s' % (sat,y,y+2),'units' : 'µm'}

            sk['%s_cot_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_cot_%s' % (sat,y)].attrs = {'long_name' : '%s cloud optical thickness averaged from %s to %s' % (sat,y,y+2),'units' : '1'}

            sk['%s_lwp_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_lwp_%s' % (sat,y)].attrs = {'long_name' : '%s liquid water path (3.7 µm) averaged from %s to %s' % (sat,y,y+2),'units' : 'g m-2'}

            sk['%s_Ctot_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_Ctot_%s' % (sat,y)].attrs = {'long_name' : '%s total cloud fraction averaged from %s to %s' % (sat,y,y+2),'units' : '1'}

            sk['%s_Clow_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_Clow_%s' % (sat,y)].attrs = {'long_name' : '%s low cloud fraction averaged from %s to %s' % (sat,y,y+2),'units' : '1'}

            sk['%s_A_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_A_%s' % (sat,y)].attrs = {'long_name' : '%s all-sky albedo averaged from %s to %s' % (sat,y,y+2),'units' : '1'}

            sk['%s_Acld_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_Acld_%s' % (sat,y)].attrs = {'long_name' : '%s overcast albedo (assuming a clear-sky albedo of 0.1) averaged from %s to %s' % (sat,y,y+2),'units' : '1'}

            sk['%s_Tskin_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_Tskin_%s' % (sat,y)].attrs = {'long_name' : '%s skin temperature averaged from %s to %s' % (sat,y,y+2),'units' : 'K'}

            sk['%s_EIS_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_EIS_%s' % (sat,y)].attrs = {'long_name' : '%s Estimated Inversion Strength averaged from %s to %s' % (sat,y,y+2),'units' : 'K'}

            sk['%s_WS_%s' % (sat,y)] = (['month','lat','lon'],np.nan*np.ones((len(sk.month),len(sk.lat),len(sk.lon))))
            sk['%s_WS_%s' % (sat,y)].attrs = {'long_name' : '%s wind speed averaged from %s to %s' % (sat,y,y+2),'units' : 'm s-1'}

            #Loop through months for temporal averaging

            ymask = np.logical_and(ds.time.dt.year>=y,ds.time.dt.year<=y+2)

            for m in range(1,13):

                mmask = ds.time.dt.month==m
                tmask = np.logical_and(ymask,mmask)

                for var1, var2 in zip(['cer','cot','lwp','Ctot','Clow','A','Acld','Tskin','EIS','WS'],['cldwatrad37_total_day_mon','cldtau_total_day_mon','lwp37_total_day_mon','C','Clow','A','Acld','aux_skint_mon','aux_inversion_mon','aux_wind_speed_mon']):
                    sk['%s_%s_%s' % (sat,var1,y)][m-1,:,:180] = np.nanmean(ds[var2][tmask,::-1,:][:,:,ds.lon>180],axis=0)
                    sk['%s_%s_%s' % (sat,var1,y)][m-1,:,180:] = np.nanmean(ds[var2][tmask,::-1,:][:,:,ds.lon<180],axis=0)
print('...Done!')


###Save file
#

filename = dir_data+'CERES/SSF1deg/SSF1deg_shipkrige_Terra.nc'
os.system('rm %s' % filename) #Delete file if it already exists
sk.to_netcdf(path=filename,mode='w')

