# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:12:24 2022

@author: DELL
"""
#%% code interpretation
# =============================================================================
# 1. This code aims to calculate partial correlation using pingouin
# 2. Compare the correlation between VPD and SIF with that of between SM and SIF
# 3. Graph spatial distribution of partial correlation in China.
# 4. Reference from https://www.machinelearningplus.com/statistics/partial-correlation/
#    https://pingouin-stats.org/generated/pingouin.partial_corr.html
# =============================================================================
#%% demo 
import pandas as pd
import matplotlib.pyplot as plt
import math
import pingouin as pg
import numpy as np
# Create a sample dataset
Data = {'A' : [4, 2, 2, 1, 8, 6, 9, 8, 11, 13, 12, 14],
        'B' : [1, 2, 2, 4, 9, 8, 9, 6, 14, 12, 13, 12],
        'Z' : [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]}        
df = pd.DataFrame(Data, columns = ['A', 'B', 'Z']) 

# Scatterplot to understand the relationship
plt.plot(df["A"],df["B"],'ro')
plt.xlabel("A")
plt.ylabel("B")

# Partial corr
pg.partial_corr(data=df, x='A', y='B', covar='Z')

# Where,
# Data = Name of the dataframe.
# x  = Name of column in dataframe.
# y = Name of column in dataframe.
# z = variable to be excluded/controlled.


# Calculate multiple variables partial corr

df.pcorr().round(7)


#%% case one for calculating partial corre between SM and VPD

## import SIF,SM,VPD datasets from google drive


## Encapsulation for graph
def Graph(lon,lat,slope1):
    from matplotlib import rcParams
    import numpy as np
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import cmaps
    import cartopy.io.shapereader as shpreader
    import matplotlib as mpl
    from cartopy.io.shapereader import Reader, natural_earth
    config = {
        "font.family": 'Times New Roman', # times new roman字体
        "font.size": 18, # 相当于小四大小
        "font.serif": ['SimSun'], # 宋体
        "mathtext.fontset": 'stix', # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
        'axes.unicode_minus': False # 处理负号，即-号
    }
    rcParams.update(config)
    plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
    
    nx, ny = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(20,20), dpi = 300)
    proj = ccrs.PlateCarree(central_longitude=80)
    # 使图一colorbar映射在同一水平
    norm1 = mpl.colors.Normalize(vmin=-abs(np.nanmax(slope1)), vmax=abs(np.nanmax(slope1)))
    ## graph one
    ax1 = fig.add_axes([0.1, 0.9, 0.5, 0.5],projection = proj)
    # ax1.set_title('(a)',loc='left')
    gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1.2, color='k', alpha=0.5, linestyle='--')
    gl.ylabels_right = False
    gl.xlabels_top = False
    gl.xlines = False
    gl.ylines = False
    china = shpreader.Reader(r'E:\Research_life\毕设分区\矢量图层-20210720T085947Z-001\矢量图层\国界与省界\bou2_4l.dbf').geometries()
    ax1.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black',zorder = 1)
    # 颜色映射，以0为分界线
    # levels = np.arange(-0.1,0.5,0.1)
    c1 = ax1.contourf(nx, ny, slope1, norm=norm1,extend = 'both',cmap = cmaps.cmp_b2r, transform = ccrs.PlateCarree())
    # c1.set_clim(vmin=-0.2, vmax=0.3)
    # c1.set_clim(vmin=-0.2, vmax=0.4)
    ax1.set_extent([73.5,135.5,18,54], crs = ccrs.PlateCarree())
    #添加南海，实际上就是新建一个子图覆盖在之前子图的右下角
    ax7 = fig.add_axes([0.5, 1.01, 0.12, 0.1], projection=ccrs.PlateCarree())
    ax7.set_extent([107, 122,3,25])
    ax7.add_geometries(Reader(r'E:\Research_life\毕设分区\矢量图层-20210720T085947Z-001\矢量图层\国界与省界\bou2_4l.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='k',linewidth=0.8)
    #添加色标，position定义色标位置，c1指定从c1填色图层取色，由于C3,C1的levles相同，所以色标一致，orientation设置色标为水平还是垂直，format设置色标标签格式
    position1 = fig.add_axes([0.1, 0.96, 0.5, 0.02])
    cbar=fig.colorbar(c1,cax=position1,orientation='horizontal',format='%.2f',)
    # cbar.set_label('$\mathrm{ΔSIF(mw/m^2/nm/sr)}$',fontsize = 24)
    return

## partial corr
SIF[np.isnan(SIF)]=0
VPD[np.isnan(VPD)]=0
SM[np.isnan(SM)]=0
NDVI_mvc[np.isnan(NDVI_mvc)]=0
VPD_yr[np.isnan(VPD_yr)]=0
SM_yr[np.isnan(SM_yr)]=0

lat=np.arange(4, 54, 0.5)    
lon=np.arange(73.5, 135.5, 0.5)


## corr between VPD and SIF
r = np.zeros((100,124))
for i in range(len(lat)):
    for j in range(len(lon)):
        try:
            df = pd.DataFrame(SIF[:,i,j],columns=['SIF'])
            df['VPD'] = pd.DataFrame(VPD[:,i,j])
            df['SM'] = pd.DataFrame(SM[:,i,j])   
            r[i,j]= pg.partial_corr(data=df, x='VPD',y='SIF',covar='SM').r
        except:
            r[i,j] = np.nan

## corr between SM and SIF
r = np.zeros((100,124))
for i in range(len(lat)):
    for j in range(len(lon)):
        try:
            df = pd.DataFrame(SIF[:,i,j],columns=['SIF'])
            df['VPD'] = pd.DataFrame(VPD[:,i,j])
            df['SM'] = pd.DataFrame(SM[:,i,j])
            r[i,j]= pg.partial_corr(data=df, x='SM',y='SIF',covar='VPD').r
        except:
            r[i,j] = np.nan


## corr between VPD and SIF
r = np.zeros((100,124))
for i in range(len(lat)):
    for j in range(len(lon)):
        df = pd.DataFrame(NDVI_mvc[:,i,j],columns=['NDVI'])
        df['VPD'] = pd.DataFrame(VPD_yr[:15,i,j])
        df['SM'] = pd.DataFrame(SM_yr[:15,i,j])   
        r[i,j]= pg.partial_corr(data=df, x='VPD',y='NDVI',covar='SM').r


## graph 
Graph(lon,lat,r)
