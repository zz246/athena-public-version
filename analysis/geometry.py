#! /usr/bin/env python2.7

def RotateEquatorToPole(lat0, lon0):
  lat1 = arcsin(cos(lat0)*cos(lon0))
  lon1 = arcsin(sin(lat0)/cos(lat1))
  if lon0 < 0. and lat0 > 0.:
    lon1 = pi - lon1
  if lon0 < 0. and lat0 < 0.:
    lon1 = - pi - lon1;
  return lat1, lon1
