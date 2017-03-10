#! /usr/bin/env python2.7
from netCDF4 import *
from pylab import *

radius = 6.371E6
data = Dataset('galewsky04-0308a.nc', 'r')

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

for t in range(ntime - 1, ntime):
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
  clines = hstack((arange(-1.1E-4, -3.E-5, 2.E-5), arange(5.E-5, 1.5E-4, 2.E-5)))
  clines = [-1.1E-4, -0.9E-4, -0.7E-4, -0.5E-4, -0.3E-4, -0.1E-4, 
            0.1E-4, 0.3E-4, 0.5E-4, 0.7E-4, 0.9E-4, 1.1E-4, 1.3E-4, 1.5E-4]
  print clines
  ax.contour(X, Y, vorticity, clines, colors = 'k', linewidths = 1.)
  #ax.contourf(X, Y, vorticity, clines)
  ax.set_ylim((0., 90.))
  #ax.set_aspect(2.)
  ax.set_aspect('equal')
  ax.set_title('144. hours')

#show()
savefig('galewsky04-0308a.png', bbox_inches = 'tight')
