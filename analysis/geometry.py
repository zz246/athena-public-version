#! /usr/bin/env python2.7
from numpy import *

def RotateEquatorToPole(lon0, lat0):
  lat1 = arcsin(cos(lat0)*cos(lon0))
  lon1 = arcsin(sin(lat0)/cos(lat1))
  try:
    if lon0 < 0. and lat0 > 0.:
      lon1 = pi - lon1
    if lon0 < 0. and lat0 < 0.:
      lon1 = - pi - lon1
  except ValueError:
    ix = where((lon0 < 0.) & (lat0 > 0.))
    lon1[ix] = pi - lon1[ix]
    ix = where((lon0 < 0.) & (lat0 < 0.))
    lon1[ix] = - pi - lon1[ix]
  return lon1, lat1

def PolarDistance(a1, r1, a2, r2):
  return sqrt((r1*cos(a1) - r2*cos(a2))**2 + (r1*sin(a1) - r2*sin(a2))**2)

def PolarCyclicDistance(angle, radius, wrap = True):
  dist = PolarDistance(angle[:-1], radius[:-1], angle[1:], radius[1:])
  last = PolarDistance(angle[-1], radius[-1], angle[0], radius[0])
  if wrap:
    return hstack((dist, last))
  else:
    return dist

def FixPolarContour(h, cutoff = 1.):
  for il, level in enumerate(h.collections):
    for kp, path in reversed(list(enumerate(level.get_paths()))):
      verts = path.vertices
      dist = PolarCyclicDistance(verts[:,0],verts[:,1])
      ix = []
      for i in range(1, len(dist)):
        if dist[i - 1] > cutoff and dist[i] > cutoff:
          ix.append(i)
      path.vertices = delete(verts, ix, axis = 0)
