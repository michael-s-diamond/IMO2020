"""
Analysis of potential cloud changes due to IMO 2020 fuel sulfur regulations

Code for making figures and calculating table values for 'Detection of large-scale cloud microphysical changes and evidence for decreasing cloud brightness within a major shipping corridor due to the 2020 International Maritime Organization marine fuel sulfur regulations'

Modifications
-------------
18 April 2023: Michael Diamond, Tallahassee, FL
    -Created
10 May 2023: Michael Diamond, Tallahassee, FL
    -Final form for ACP Letters initial submission
"""

#Import libraries
import numpy as np
import xarray as xr
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.special import expit
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from glob import glob
import os
import warnings

#Set paths
dir_data = '/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/Results/'
dir_SSF = '/Users/michaeldiamond/Documents/Data/CERES/SSF1deg/'
dir_figs = '/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/Figures/'

#Load data
dKr = {} #Kriging results
for var in ['cer','Acld']:
    dKr[var] = {}
    for y in list(range(2002,2023,3))+['clim']:
        dKr[var][str(y)] = {}
        for sea, M in zip(['SON','ANN'],['M9to11','M1to12']):
            dKr[var][str(y)][sea] = xr.open_dataset(dir_data+'Data_Terra_%s_%s_C_%s.nc' % (var,y,M))
Terra = xr.open_dataset(dir_SSF+'SSF1deg_shipkrige_Terra.nc') #Kriging input
SSF = xr.open_dataset(glob(dir_SSF+'CERES_SSF1deg-Month_Terra-MODIS_Ed4.1_Subset_*.nc')[0]) #SSF1deg data


#Extract useful information
trk = np.zeros(70,dtype=bool)
trk[2::7] = True #
trk[3::7] = True ##Central track only
trk[4::7] = True #

dEst, dObs = {}, {} #Estimated and observed central track values

for var in ['cer','Acld']:
    dEst[var], dObs[var] = {}, {}
    for y in list(range(2002,2023,3))+['clim']:
        dEst[var][str(y)], dObs[var][str(y)] = {}, {}
        for sea in ['SON','ANN']:
            if var == 'cer':
                dEst[var][str(y)][sea] = dKr[var][str(y)][sea]['krSims'][:,trk].mean(axis=-1)
                dObs[var][str(y)][sea] = dKr[var][str(y)][sea]['shipObs'][trk].mean()
            elif var == 'Acld':
                dEst[var][str(y)][sea] = expit(dKr[var][str(y)][sea]['krSims'][:,trk]).mean(axis=-1)
                dObs[var][str(y)][sea] = expit(dKr[var][str(y)][sea]['shipObs'][trk]).mean()
                

#Significance testing
praw = []
for var in ['cer','Acld']:
    for y in list(range(2002,2023,3))+['clim']:
        for sea in ['SON','ANN']:
            praw.append(dKr[var][str(y)][sea].pVal)
padj = multipletests(praw,alpha=.05,method='fdr_bh')[1] #BH correction for mulitple testing
i = 0
for var in ['cer','Acld']:
    for y in list(range(2002,2023,3))+['clim']:
        for sea in ['SON','ANN']:
            dKr[var][str(y)][sea].attrs['padj'] = padj[i]
            i += 1

"""
Figures 1, 2, S1, S2: Maps of SON cer, SON Acld, ANN cer, and ANN Acld values
"""

def map_plot(var='cer',sea='SON'):
    """
    Plot Ship, NoShip, and Difference values on a map projection centered on the SEA
    
    Parameters
    ----------
    var : str
    'cer' for effective radius and 'Acld' for overcast albedo
    
    sea : str
    'SON' for austral spring and 'ANN' for annual average
    """

    #Set up plot
    fig = plt.figure(figsize=(13,9))
    plt.clf()
    
    fs = 16
    plabs = ['a','b','c','d','e','f','g','h','i']
    
    #Plot data on 3x3 grid of map panels
    n = 1
    for y in ['clim','2017','2020']:
        for kind in ['Ship','NoShip','Difference']:

            #Get data and set value limits, colors, etc.
            if kind == 'Ship':
                if var == 'cer':
                    data = (dKr[var][str(y)][sea]['Obs']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    cmap = 'viridis_r'
                    vmin = 10
                    vmax = 16
                    ticks = np.arange(10,16.1,2)
                    lab = r'$r_\mathrm{e,Ship}$ (µm)'
                elif var == 'Acld':
                    data = 100*expit(dKr[var][str(y)][sea]['Obs']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    cmap = 'viridis'
                    if sea == 'SON':
                        vmin = 28
                        vmax = 38
                        ticks = np.arange(28,38.1,2)
                    elif sea == 'ANN':
                        vmin = 24
                        vmax = 34
                        ticks = np.arange(24,34.1,2)
                    lab = r'$A_\mathrm{cld,Ship}$ (%)'
            if kind == 'NoShip':
                if var == 'cer':
                    data = (dKr[var][str(y)][sea]['Est']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    cmap = 'viridis_r'
                    vmin = 10
                    vmax = 16
                    ticks = np.arange(10,16.1,2)
                    ticks = np.arange(10,16.1,2)
                    lab = r'$r_\mathrm{e,NoShip}$ (µm)'
                elif var == 'Acld':
                    data = 100*expit(dKr[var][str(y)][sea]['Est']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    cmap = 'viridis'
                    if sea == 'SON':
                        vmin = 28
                        vmax = 38
                        ticks = np.arange(28,38.1,2)
                    elif sea == 'ANN':
                        vmin = 24
                        vmax = 34
                        ticks = np.arange(24,34.1,2)
                    lab = r'$A_\mathrm{cld,Ship}$ (%)'
            if kind == 'Difference':
                if var == 'cer':
                    data = (dKr[var][str(y)][sea]['Obs']-dKr[var][str(y)][sea]['Est']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    cmap = 'RdYlBu_r'
                    if sea == 'SON':
                        vmin = -.8
                        vmax = .8
                        ticks = np.arange(-.8,.81,.4)
                    elif sea == 'ANN':
                        vmin = -.5
                        vmax = .5
                        ticks = np.arange(-.5,.51,.25)
                    lab = r'$\Delta r_\mathrm{e}$ (µm)'
                    high = dKr[var][str(y)][sea]['highEst'].sel(lat=slice(-25,-2),lon=slice(-15,15))
                    low = dKr[var][str(y)][sea]['lowEst'].sel(lat=slice(-25,-2),lon=slice(-15,15))
                    obs = dKr[var][str(y)][sea]['Obs'].sel(lat=slice(-25,-2),lon=slice(-15,15))
                    sig = np.logical_or(obs>high,obs<low)
                    pval = dKr[var][str(y)][sea].padj
                elif var == 'Acld':
                    data = 100*(expit(dKr[var][str(y)][sea]['Obs'])-expit(dKr[var][str(y)][sea]['Est'])).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    cmap = 'RdYlBu'
                    if sea == 'SON':
                        vmin = -2
                        vmax = 2
                        ticks = np.arange(-2,2.1,1)
                    elif sea == 'ANN':
                        vmin = -1
                        vmax = 1
                        ticks = np.arange(-1,1.1,.5)
                    lab = r'$\Delta A_\mathrm{cld}$ (%)'
                    high = expit(dKr[var][str(y)][sea]['highEst']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    low = expit(dKr[var][str(y)][sea]['lowEst']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    obs = expit(dKr[var][str(y)][sea]['Obs']).sel(lat=slice(-25,-2),lon=slice(-15,15))
                    sig = np.logical_or(obs>high,obs<low)
                    pval = dKr[var][str(y)][sea].padj

            #Set up map
            lon, lat = np.meshgrid(data.lon,data.lat)

            subplt = 330+n

            ax = plt.subplot(subplt, projection=ccrs.Mercator())
            ax.set_extent([-15,15,-19.5,-6])
            ax.add_feature(cartopy.feature.LAND,facecolor='k',edgecolor='w',lw=.5,zorder=10)
            ax.add_feature(cartopy.feature.BORDERS,edgecolor='w',lw=.5,zorder=10)

            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='.5', alpha=0,zorder=9)
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xlocator = mticker.FixedLocator([-13,0,8])
            gl.ylocator = mticker.FixedLocator([-18,-13,-8])
            gl.xlabel_style = {'size': fs-4}
            gl.ylabel_style = {'size': fs-4}
            gl.yformatter = LatitudeFormatter(degree_symbol='° ')
            gl.xformatter = LongitudeFormatter(degree_symbol='° ')

            if kind == 'Difference': #Significance test results
                ax.scatter(lon[sig.values],lat[sig.values],transform=ccrs.PlateCarree(),c='w',s=2,zorder=11)
                if pval < .001:
                    ax.text(-12.5,-17.5,r'$p_\mathrm{field} \ll$0.001',fontsize=10,ha='left',va='bottom',transform=ccrs.PlateCarree())
                else:
                    ax.text(-12.5,-17.5,r'$p_\mathrm{field}$=%.3f' % pval,fontsize=10,ha='left',va='bottom',transform=ccrs.PlateCarree())

            #Plot values
            shading = ax.pcolormesh(lon,lat,data.values,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)

            cbar = plt.colorbar(shading,extend='both',orientation='horizontal',ticks=ticks)
            if n > 6: cbar.set_label(lab,fontsize=fs)
            cbar.ax.tick_params(labelsize=fs-2)

            ax.plot([-13,8,8,-13,-13],[-18,-18,-8,-8,-18],'k',lw=1,transform=ccrs.PlateCarree())

            #Panel labels
            ax.text(14.95,-6.325,'(%s)' % plabs[n-1],color='w',fontsize=10,transform=ccrs.PlateCarree(),zorder=11,va='top',ha='right') 

            #Row/column labels
            if n == 1: 
                plt.title('Ship\n',fontsize=fs)
                ax.text(-.3,.5,s='2002–2019',transform = ax.transAxes,fontsize=fs,va='center',ha='center',rotation=90)

            elif n == 2: plt.title('NoShip\n',fontsize=fs)
            elif n == 3: plt.title('Difference\n',fontsize=fs)
            elif n == 4: ax.text(-.3,.5,s='2017–2019',transform = ax.transAxes,fontsize=fs,va='center',ha='center',rotation=90)
            elif n == 7: ax.text(-.3,.5,s='2020–2022',transform = ax.transAxes,fontsize=fs,va='center',ha='center',rotation=90)

            n += 1

    plt.savefig(dir_figs+'maps_%s_%s_clim_2017_2020.png' % (var,sea),dpi=450)


"""
Figures S3 and S4: Maps of SON and ANN meteorology
"""

SSTv_SON = np.zeros((6,23,30))
SSTv_ANN = np.zeros((6,23,30))
EISv_SON = np.zeros((6,23,30))
EISv_ANN = np.zeros((6,23,30))
WSv_SON = np.zeros((6,23,30))
WSv_ANN = np.zeros((6,23,30))

for i, y in zip(range(6),range(2002,2020,3)):
    SSTv_SON[i][:] = Terra['Terra_Tskin_%s' % y].sel(month=slice(9,11),lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0).values
    SSTv_ANN[i][:] = Terra['Terra_Tskin_%s' % y].sel(lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0).values
    EISv_SON[i][:] = Terra['Terra_EIS_%s' % y].sel(month=slice(9,11),lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0).values
    EISv_ANN[i][:] = Terra['Terra_EIS_%s' % y].sel(lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0).values
    WSv_SON[i][:] = Terra['Terra_WS_%s' % y].sel(month=slice(9,11),lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0).values
    WSv_ANN[i][:] = Terra['Terra_WS_%s' % y].sel(lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0).values

dVar = {'SON' : {}, 'ANN' : {}}
dVar['SON']['Tskin'] = np.std(SSTv_SON,axis=0)
dVar['ANN']['Tskin'] = np.std(SSTv_ANN,axis=0)
dVar['SON']['EIS'] = np.std(EISv_SON,axis=0)
dVar['ANN']['EIS'] = np.std(EISv_ANN,axis=0)
dVar['SON']['WS'] = np.std(WSv_SON,axis=0)
dVar['ANN']['WS'] = np.std(WSv_ANN,axis=0)
    
def met_plot(sea='SON'):
    with warnings.catch_warnings(): #Ignore runtime warnings
        warnings.simplefilter("ignore")
    
        #Set up plot
        fig = plt.figure(figsize=(13,9))
        plt.clf()

        fs = 16
        plabs = ['a','b','c','d','e','f','g','h','i']

        #Plot data on 3x3 grid of map panels
        n = 1
        for y in ['clim','2017','2020']:
            for var in ['Tskin','EIS','WS']:

                if y == 'clim': 
                    if sea == 'SON':
                        data = Terra['Terra_'+var+'_'+y].sel(month=slice(9,11),lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0)
                    elif sea == 'ANN':
                        data = Terra['Terra_'+var+'_'+y].sel(lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0)
                    cmap = 'viridis'
                    if var == 'Tskin':
                        if sea == 'SON':
                            vmin = 290
                            vmax = 300
                        elif sea == 'ANN':
                            vmin = 290
                            vmax = 300
                        ticks = np.arange(vmin,vmax+.1,2)
                    elif  var == 'EIS':
                        if sea == 'SON':
                            vmin = 2
                            vmax = 12
                        elif sea == 'ANN':
                            vmin = 0
                            vmax = 10
                        ticks = np.arange(vmin,vmax+.1,2)
                    elif var == 'WS':
                        if sea == 'SON':
                            vmin = 3
                            vmax = 8
                        elif sea == 'ANN':
                            vmin = 2
                            vmax = 7
                        ticks = np.arange(vmin,vmax+.1,1)

                else: 
                    if sea == 'SON': 
                        data = Terra['Terra_'+var+'_'+y].sel(month=slice(9,11),lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0)-Terra['Terra_'+var+'_clim'].sel(month=slice(9,11),lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0)
                                            
                    elif sea == 'ANN': 
                        data = Terra['Terra_'+var+'_'+y].sel(lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0)-Terra['Terra_'+var+'_clim'].sel(lat=slice(-2,-25),lon=slice(-15,15)).mean(axis=0)
                    cmap = 'RdYlBu_r'

                    if var == 'Tskin':
                        vmin = -1
                        vmax = 1
                        ticks = np.arange(vmin,vmax+.1,.5)
                        lab = r'$T_\mathrm{skin}$ (K)'
                    elif  var == 'EIS':
                        vmin = -1
                        vmax = 1
                        ticks = np.arange(vmin,vmax+.1,.5)
                        lab = r'EIS (K)'
                    elif var == 'WS':
                        vmin = -.5
                        vmax = .5
                        ticks = np.arange(vmin,vmax+.1,.25)
                        lab = r'Wind speed ($\mathrm{m}$ $\mathrm{s}^{-1}$)'


                #Set up map
                lon, lat = np.meshgrid(data.lon,data.lat)

                subplt = 330+n

                ax = plt.subplot(subplt, projection=ccrs.Mercator())
                ax.set_extent([-15,15,-19.5,-6])
                ax.add_feature(cartopy.feature.LAND,facecolor='k',edgecolor='w',lw=.5,zorder=10)
                ax.add_feature(cartopy.feature.BORDERS,edgecolor='w',lw=.5,zorder=10)

                gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=1, color='.5', alpha=0,zorder=9)
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xlocator = mticker.FixedLocator([-13,0,8])
                gl.ylocator = mticker.FixedLocator([-18,-13,-8])
                gl.xlabel_style = {'size': fs-4}
                gl.ylabel_style = {'size': fs-4}
                gl.yformatter = LatitudeFormatter(degree_symbol='° ')
                gl.xformatter = LongitudeFormatter(degree_symbol='° ')

                #Plot values
                shading = ax.pcolormesh(lon,lat,data.values,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
                
                #Stippling for anomalies > 2 sigma
                if y != 'clim':
                    dvar = dVar[sea][var]
                    sig = np.abs(data) > 2*dvar
                    ax.scatter(lon[sig.values],lat[sig.values],transform=ccrs.PlateCarree(),c='w',s=2,zorder=9)

                cbar = plt.colorbar(shading,extend='both',orientation='horizontal',ticks=ticks)
                if n > 6: cbar.set_label(lab,fontsize=fs)
                cbar.ax.tick_params(labelsize=fs-2)

                ax.plot([-13,8,8,-13,-13],[-18,-18,-8,-8,-18],'k',lw=1,transform=ccrs.PlateCarree())

                #Panel labels
                ax.text(14.95,-6.325,'(%s)' % plabs[n-1],color='w',fontsize=10,transform=ccrs.PlateCarree(),zorder=11,va='top',ha='right') 

                #Row/column labels
                if n == 1: 
                    plt.title('SST\n',fontsize=fs)
                    ax.text(-.3,.5,s='2002–2019\nmean',transform = ax.transAxes,fontsize=fs,va='center',ha='center',rotation=90)

                elif n == 2: plt.title('EIS\n',fontsize=fs)
                elif n == 3: plt.title('WS\n',fontsize=fs)
                elif n == 4: ax.text(-.3,.5,s='2017–2019\nanomaly',transform = ax.transAxes,fontsize=fs,va='center',ha='center',rotation=90)
                elif n == 7: ax.text(-.3,.5,s='2020–2022\nanomaly',transform = ax.transAxes,fontsize=fs,va='center',ha='center',rotation=90)

                n += 1

            plt.savefig(dir_figs+'maps_met_%s_clim_2017_2020.png' % sea,dpi=450)
    

"""
Table S1: Summary statistics
"""

#Mean ship values

def print_row(var='cer',y='clim',sea='SON'):
    print(var,sea,y)
    
    if var == 'cer':
        ship_ = dObs[var][y][sea].values
        noship_ = dEst[var][y][sea].values
    elif var == 'Acld':
        ship_ = 100*dObs[var][y][sea].values
        noship_ = 100*dEst[var][y][sea].values

    abs_ = np.mean(ship_-noship_)
    abs_l = np.percentile(ship_-noship_,2.5)
    abs_h = np.percentile(ship_-noship_,97.5)
        
    rel_ = np.mean(100*(ship_-noship_)/ship_)
    rel_l = np.percentile(100*(ship_-noship_)/ship_,2.5)
    rel_h = np.percentile(100*(ship_-noship_)/ship_,97.5)
        
    p_ = dKr[var][y][sea].padj
        
    reg_ = dKr[var][y][sea].parSel
     
    print('Mean ship:','%.2f' % ship_)
    print('Abs diff:','%.2f (%.2f to %.2f)' % (abs_,abs_l,abs_h))
    print('Rel diff:','%.1f (%.1f to %.1f)' % (rel_,rel_l,rel_h))
    print('p:','%.5f' % p_)
    print('regs:',reg_)


"""
Figure 3: PDFs
"""

plt.figure(figsize=(9,14))
plt.clf()
fs = 16

colors = {'clim' : '.5', '2002' : cm.Spectral(.05), '2005' : cm.Spectral(.2), '2008' : cm.Spectral(.33), '2011' : cm.Spectral(.67), '2014' : cm.Spectral(.8), '2017' : cm.Spectral(.95), '2020' : 'k'}
titles = [r'(a) Relative change in austral spring $r_\mathrm{e}$ (%)',r'(b) Relative change in austral spring $A_\mathrm{cld}$ (%)',r'(c) Relative change in annual mean $r_\mathrm{e}$ (%)',r'(d) Relative change in annual mean $A_\mathrm{cld}$ (%)']

n = 1
for var, sea, xvals in zip(2*['cer','Acld'],2*['SON']+2*['ANN'],[np.arange(-8.25,2.75,.01),np.arange(-1.25,5.75,.01),np.arange(-4.75,.25,.01),np.arange(-.75,4.25,.01)]):

    plt.subplot(410+n)

    ymax = 0

    for y in list(range(2002,2023,3))+['clim']:

        y = str(y)

        data = 100*(dObs[var][y][sea].values-dEst[var][y][sea].values)/dObs[var][y][sea].values

        gdensity = gaussian_kde(data)
        yline = gdensity(xvals)
        if ymax < np.max(yline): ymax = np.max(yline)

        pval = dKr[var][y][sea].padj
        if pval < .0001:
            ls = 'solid'
        elif pval > .0001 and pval < .01:
            ls = 'dashed'
        elif pval > .01 and pval < .1:
            ls = 'dotted'

        if y == 'clim':
            plt.plot(xvals,yline,c=colors[y],lw=3,zorder=1,ls=ls)
            plt.fill_between(xvals,0*xvals,yline,facecolor='.75')
        elif y == '2020':
            plt.plot(xvals,yline,c=colors[y],lw=3,zorder=10,ls=ls)
        else:
            plt.plot(xvals,yline,c=colors[y],lw=3,zorder=5,ls=ls)

        if var == 'cer' and sea == 'SON':
            if y == 'clim': plt.fill_between(np.nan*xvals,np.nan*xvals,np.nan*xvals,facecolor='.75',label='2002–2019')
            else: plt.plot(np.nan*xvals,np.nan*xvals,c=colors[y],label='%s–%s' % (y,int(y)+2))

    if var == 'Acld' and sea == 'SON':
        plt.plot(np.nan*xvals,np.nan*xvals,c='k',label=r'$p_\mathrm{field}$<0.0001',ls='solid')
        plt.plot(np.nan*xvals,np.nan*xvals,c='k',label=r'0.0001<$p_\mathrm{field}$<0.01',ls='dashed')
        plt.plot(np.nan*xvals,np.nan*xvals,c='k',label=r'0.01<$p_\mathrm{field}$<0.1',ls='dotted')


    plt.plot(xvals,0*xvals,'k',lw=2,zorder=12)
    plt.plot([0,0],[0,2*ymax],'k--',lw=1)

    if sea == 'SON':
        plt.legend(frameon=False,fontsize=fs-2,loc=1)

    plt.ylim(0,1.05*ymax)
    plt.xlim(np.min(xvals),np.max(xvals))

    plt.yticks(fontsize=fs-2)
    plt.ylabel('Probability density',fontsize=fs)
    plt.xticks(fontsize=fs-2)

    plt.title(titles[n-1],loc='left',fontsize=fs)

    n += 1

plt.tight_layout()

plt.savefig(dir_figs+'PDFs.png',dpi=450)
plt.savefig(dir_figs+'PDFs.eps',dpi=450)


"""
Table S2: Percentiles of 2020/clim ratios
"""

def print_ratio(var='cer',sea='SON'):
    print(var,sea)
    
    d20 = (dObs[var]['2020'][sea]-dEst[var]['2020'][sea])/dObs[var]['2020'][sea]
    dcl = (dObs[var]['clim'][sea]-dEst[var]['clim'][sea])/dObs[var]['clim'][sea]
    
    ratio = d20/dcl
    
    print('%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f' % (np.percentile(ratio,0.5),np.percentile(ratio,2.5),np.percentile(ratio,5),np.percentile(ratio,25),np.percentile(ratio,50),np.percentile(ratio,75),np.percentile(ratio,95),np.percentile(ratio,97.5),np.percentile(ratio,99.5)))


"""
Figure 4: Time series
"""

plt.figure(figsize=(10,7.5))
plt.clf()
fs = 16

var = 'cer'
sea = 'SON'

tsla_Obs = np.array([dKr[var][y][sea].sel(lat=slice(-18,-8),lon=slice(-13,8))['Obs'].mean().values for y in list(dObs[var].keys())[:-1]])
tsla_Est = np.array([dKr[var][y][sea].sel(lat=slice(-18,-8),lon=slice(-13,8))['Est'].mean().values for y in list(dObs[var].keys())[:-1]])

tscl_Obs = np.array([dKr[var][y][sea].sel(lat=slice(-18,-8),lon=slice(-13,8))['Obs'].mean().values for y in list(dObs[var].keys())[-1:]])
tscl_Est = np.array([dKr[var][y][sea].sel(lat=slice(-18,-8),lon=slice(-13,8))['Est'].mean().values for y in list(dObs[var].keys())[-1:]])

alt2020 = tsla_Est[-1]+(tscl_Obs-tscl_Est)

plt.scatter(np.arange(2003.5,2022,3),tsla_Obs,c='k',label='Ship')
plt.plot([2002,2005,2005,2008,2008,2011,2011,2014,2014,2017,2017,2020,2020,2023],np.array([(tsla_Obs[i],tsla_Obs[i]) for i in range(len(tsla_Obs))]).ravel(),c='k',lw=1)

plt.scatter(np.arange(2003.5,2022,3),tsla_Est,color=cm.Blues(.5),marker='D',label='NoShip')
plt.plot([2002,2005,2005,2008,2008,2011,2011,2014,2014,2017,2017,2020,2020,2023],np.array([(tsla_Est[i],tsla_Est[i]) for i in range(len(tsla_Obs))]).ravel(),c=cm.Blues(.5),lw=1,ls='dashed')

plt.scatter(2021.5,alt2020,c='k',marker='x',label='Noncompliance')

plt.plot(2*[2021.5+.25],[alt2020,tsla_Obs[-1]],c=cm.Reds(.75),ls='dotted',lw=2,label='True IMO effect,\nkriging method')
plt.text(2021.5+.25,(tsla_Obs[-1]+alt2020)/2,s='}',fontsize=24,ha='left',va='center',color=cm.Reds(.75),)
plt.text(2021.5+.65,(tsla_Obs[-1]+alt2020)/2,s='%.1f\nµm' % (tsla_Obs[-1]-alt2020),fontsize=12,ha='left',va='center',color=cm.Reds(.75))

plt.plot(2*[2021.5-.25],[tsla_Obs[-2],alt2020],c=cm.Reds(.5),ls='dotted',lw=2,label='False IMO effect,\npersistence method')
plt.text(2021.5-.25,(alt2020+tsla_Obs[-2])/2,s='{',fontsize=24,ha='right',va='center',color=cm.Reds(.5))
plt.text(2021.5-.65,(alt2020+tsla_Obs[-2])/2,s='%.1f\nµm' % (alt2020-tsla_Obs[-2]),fontsize=12,ha='right',va='center',color=cm.Reds(.5))

plt.ylabel(r'Austral spring regional mean $r_\mathrm{e}$ (µm)',fontsize=fs)
plt.yticks(fontsize=fs-2)
plt.xticks(np.arange(2002,2024,3),fontsize=fs-2)

plt.xlim(2002,2023)
plt.ylim(12.5,13.7)

plt.legend(frameon=False,fontsize=fs-4)

plt.savefig(dir_figs+'NoIMO.png',dpi=450)
plt.savefig(dir_figs+'NoIMO.eps',dpi=450)


"""
Figure 5: IRFaci calculations
"""

dIRF = {'cer' : {'clim' : {}, '2020' : {}}, 'Acld' : {'clim' : {}, '2020' : {}}}

phi = 0.6
acld = 0.5

#
###SON clim
#
sea = 'SON'
yy = 'clim'

S = SSF['toa_solar_all_mon'][np.logical_and(np.logical_and(SSF.time.dt.year>=2002,SSF.time.dt.year<=2019),np.logical_and(SSF.time.dt.month>=9,SSF.time.dt.month<=11))].sel(lat=slice(-18,-8))[:,:,np.logical_or(SSF.lon<=8,SSF.lon>=360-13)].mean().values
C = Terra['Terra_Clow_'+yy].sel(month=slice(9,11),lat=slice(-8,-18),lon=slice(-13,8)).mean().values

dIRF['cer'][yy][sea] = -S*C*phi*acld*(1-acld)*-(dObs['cer'][yy][sea]-dEst['cer'][yy][sea])/dObs['cer'][yy][sea]
dIRF['Acld'][yy][sea] = -S*C*(dObs['Acld'][yy][sea]-dEst['Acld'][yy][sea])

#
###SON 2020
#
sea = 'SON'
yy = '2020'

S = SSF['toa_solar_all_mon'][np.logical_and(np.logical_and(SSF.time.dt.year>=2020,SSF.time.dt.year<=2022),np.logical_and(SSF.time.dt.month>=9,SSF.time.dt.month<=11))].sel(lat=slice(-18,-8))[:,:,np.logical_or(SSF.lon<=8,SSF.lon>=360-13)].mean().values
C = Terra['Terra_Clow_'+yy].sel(month=slice(9,11),lat=slice(-8,-18),lon=slice(-13,8)).mean().values

dIRF['cer'][yy][sea] = -S*C*phi*acld*(1-acld)*-(dObs['cer'][yy][sea]-dEst['cer'][yy][sea])/dObs['cer'][yy][sea]
dIRF['Acld'][yy][sea] = -S*C*(dObs['Acld'][yy][sea]-dEst['Acld'][yy][sea])

#
###ANN clim
#
sea = 'ANN'
yy = 'clim'

S = SSF['toa_solar_all_mon'][np.logical_and(SSF.time.dt.year>=2002,SSF.time.dt.year<=2019)].sel(lat=slice(-18,-8))[:,:,np.logical_or(SSF.lon<=8,SSF.lon>=360-13)].mean().values
C = Terra['Terra_Clow_'+yy].sel(lat=slice(-8,-18),lon=slice(-13,8)).mean().values

dIRF['cer'][yy][sea] = -S*C*phi*acld*(1-acld)*-(dObs['cer'][yy][sea]-dEst['cer'][yy][sea])/dObs['cer'][yy][sea]
dIRF['Acld'][yy][sea] = -S*C*(dObs['Acld'][yy][sea]-dEst['Acld'][yy][sea])

#
###ANN 2020
#
sea = 'ANN'
yy = '2020'

S = SSF['toa_solar_all_mon'][np.logical_and(SSF.time.dt.year>=2020,SSF.time.dt.year<=2022)].sel(lat=slice(-18,-8))[:,:,np.logical_or(SSF.lon<=8,SSF.lon>=360-13)].mean().values
C = Terra['Terra_Clow_'+yy].sel(lat=slice(-8,-18),lon=slice(-13,8)).mean().values

dIRF['cer'][yy][sea] = -S*C*phi*acld*(1-acld)*-(dObs['cer'][yy][sea]-dEst['cer'][yy][sea])/dObs['cer'][yy][sea]
dIRF['Acld'][yy][sea] = -S*C*(dObs['Acld'][yy][sea]-dEst['Acld'][yy][sea])


#
###Plot
#
plt.figure(figsize=(9,7.5))
plt.clf()
fs = 16

xvals = np.arange(-5,5,.01)


plt.subplot(2,1,1) #SON

plt.fill_between(xvals,0*xvals,gaussian_kde(dIRF['cer']['clim']['SON'])(xvals),facecolor=cm.Blues(.1),label='2002–2019')
plt.fill_between(xvals,0*xvals,gaussian_kde(dIRF['cer']['2020']['SON'])(xvals),facecolor='.8',label='2020–2022')
plt.fill_between(xvals,0*xvals,gaussian_kde(dIRF['cer']['2020']['SON']-dIRF['cer']['clim']['SON'])(xvals),facecolor=cm.Reds(.3),label='IMO 2020')

plt.plot(xvals,gaussian_kde(dIRF['Acld']['clim']['SON'])(xvals),c=cm.Blues(.5),lw=2)
plt.plot(xvals,gaussian_kde(dIRF['Acld']['2020']['SON'])(xvals),c='k',lw=2)
plt.plot(xvals,gaussian_kde(dIRF['Acld']['2020']['SON']-dIRF['Acld']['clim']['SON'])(xvals),c=cm.Reds(.75),lw=2)

plt.plot([0,0],[-10,10],'k--',lw=1)

plt.legend(frameon=False,fontsize=fs-2,loc=1)

plt.ylim(0,2.75)
plt.xlim(-4.75,4.5)

plt.xticks(np.arange(-4,5,1),fontsize=fs-2)
plt.ylabel('Probability density',fontsize=fs)
plt.yticks(fontsize=fs-2)

plt.title(r'(a) Austral spring Twomey effect estimates ($\mathrm{W}$ $\mathrm{m}^{-2}$)',fontsize=fs,loc='left')


plt.subplot(2,1,2) #ANN

plt.fill_between(xvals,0*xvals,gaussian_kde(dIRF['cer']['clim']['ANN'])(xvals),facecolor=cm.Blues(.1))
plt.fill_between(xvals,0*xvals,gaussian_kde(dIRF['cer']['2020']['ANN'])(xvals),facecolor='.8',label=r'$-{F}_{\odot} \mathcal{C}_\mathrm{low} \phi_{\mathrm{atm}} \alpha_\mathrm{cld} (1-\alpha_\mathrm{cld}) \frac{-\Delta {r}_\mathrm{e}}{r_\mathrm{e,Ship}}$')
plt.fill_between(xvals,0*xvals,gaussian_kde(dIRF['cer']['2020']['ANN']-dIRF['cer']['clim']['ANN'])(xvals),facecolor=cm.Reds(.3))

plt.plot(xvals,gaussian_kde(dIRF['Acld']['clim']['ANN'])(xvals),c=cm.Blues(.5),lw=2)
plt.plot(xvals,gaussian_kde(dIRF['Acld']['2020']['ANN'])(xvals),c='k',lw=2,label=r'$-{F}_{\odot} \mathcal{C}_\mathrm{low} \Delta {A}_\mathrm{cld}$')
plt.plot(xvals,gaussian_kde(dIRF['Acld']['2020']['ANN']-dIRF['Acld']['clim']['ANN'])(xvals),c=cm.Reds(.75),lw=2)

plt.plot([0,0],[-10,10],'k--',lw=1)

plt.legend(frameon=False,fontsize=fs-2,loc=1)

plt.ylim(0,3.33)
plt.xlim(-4.75,4.5)

plt.xticks(np.arange(-4,5,1),fontsize=fs-2)
plt.ylabel('Probability density',fontsize=fs)
plt.yticks(fontsize=fs-2)

plt.title(r'(b) Annual mean Twomey effect estimates ($\mathrm{W}$ $\mathrm{m}^{-2}$)',fontsize=fs,loc='left')

plt.tight_layout()

plt.savefig(dir_figs+'Twomey.png',dpi=450)
plt.savefig(dir_figs+'Twomey.eps')


"""
Figures S5-8: Semivariograms
"""

h = np.linspace(0,22,1000)
def exp_semivar(x,sigma2,phi):
    return sigma2*(1-np.exp(-h/phi))


def plot_semivar(var='cer',sea='SON'):

    plabs = ['a','b','c','d','e','f','g','h']
    
    plt.figure(figsize=(6.5,8))
    plt.clf()
    fs = 10

    n = 1
    
    for y in list(dObs[var].keys()):

        plt.subplot(4,2,n)

        sigma2 = dKr[var][y][sea].Sigma2
        phi = dKr[var][y][sea].Phi
        bins = dKr[var][y][sea]['Semivariance']

        plt.plot(h,exp_semivar(h,sigma2,phi),c='c',lw=2,label='Fitted semivariogram')
        plt.scatter(bins.bins,bins.values,color='k',s=25,label='Binned semivariogram')

        plt.xlim(-.5,22)
        plt.ylim(0,1.1*np.max(bins.values))

        plt.yticks(np.linspace(0,np.max(bins),3),['%.1e' % yt for yt in np.linspace(0,np.max(bins),3)],fontsize=fs-2)
        plt.xticks(fontsize=fs-2)
        
        if n%2 == 1: plt.ylabel('Semivariance',fontsize=fs)
        if n > 6: plt.xlabel('Distance (grid boxes)',fontsize=fs)
        if n == 8: 
            if sea == 'ANN': plt.legend(frameon=False,fontsize=fs-2,loc=2)
            else: plt.legend(frameon=False,fontsize=fs-2,loc=4)

        if var == 'cer':
            name = r'$r_\mathrm{e}$ (µm)'
        elif var == 'Acld':
            name = r'logit($A_\mathrm{cld}$)'

        if y == 'clim':
            time = '2002–2019'
        else: time = '%i–%i' % (int(y),int(y)+2)

        plt.title('(%s) %s %s %s' % (plabs[n-1], time, sea, name),fontsize=fs,loc='left')
        
        n += 1
    
    plt.tight_layout()

    plt.savefig(dir_figs+'semivar_%s_%s_%s.png' % (var,sea,y),dpi=450)
    plt.savefig(dir_figs+'semivar_%s_%s_%s.eps' % (var,sea,y))

