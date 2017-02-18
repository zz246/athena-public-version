#! /usr/bin/env python2.7
from netCDF4 import *
from pylab import *

radius = 6.371E6
data = Dataset('galewsky04-0214c.nc', 'r')

time = data.variables['time'][:]
lon = data.variables['x1'][:]
lat = data.variables['x2'][:]

ntime, nlat, nlon = len(time), len(lat), len(lon)
print ntime, nlat, nlon

hgt = squeeze(data.variables['rho'][:])
uwind = squeeze(data.variables['vel1'][:])
vwind = squeeze(data.variables['vel2'][:])

data.close()

vorticity = zeros((nlat - 2, nlon - 2))

X, Y = meshgrid(lon[1:-1] * 180./pi, lat[1:-1] * 180./pi)

for t in range(433, 434):
  print time[t]/3600.
#for t in range(1):
  figure(1, figsize = (12, 8))
  ax = axes()
  for j in range(1, nlat - 1):
    for i in range(1, nlon - 1):
      vorticity[j - 1, i - 1] =  1./(radius*cos(lat[j])) * ( \
        (vwind[t, j, i + 1] - vwind[t, j, i - 1]) / (lon[i + 1] - lon[i - 1]) \
        - (cos(lat[j+1])*uwind[t, j + 1, i] - cos(lat[j-1])*uwind[t, j - 1, i]) / (lat[j + 1] - lat[j - 1]))
  print vorticity.max()
  print vorticity.min()
  clines = hstack((arange(-16.E-5, 0, 2.E-5), arange(2.E-5, 16.E-5, 2.E-5)))
  ax.contour(X, Y, vorticity, clines, colors = 'k', linewidths = 1.)
  ax.set_ylim((10., 80.))
  ax.set_aspect(2.)
  ax.set_title('144. hours')

savefig('galewsky04-0214c.png', bbox_inches = 'tight')

