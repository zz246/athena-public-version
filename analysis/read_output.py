#! /usr/bin/env python2.7
from glob import glob
from numpy import genfromtxt, zeros, vstack 
import re

def GetBlockId(fname):
  for name in fname.split('.'):
    if name[:5] == 'block':
      return int(name[5:])
  raise Exception('Block id not found')

def GetFrameId(fname):
  fid = fname.split('.')[-2]
  return int(fid)

def GetCasename(fname):
  return fname.split('.')[0]

def GetExtension(fname):
  return fname.split('.')[-1]

def GetOutputId(fname):
  for name in fname.split('.'):
    if name[:3] == 'out':
      return int (name[3:])
  raise Exception('Output id not found')

def GetNumberBlocks(ext, path = './'):
  fnames = glob('*.' + ext)
  return len(set(map(GetBlockId, fnames)))

def GetNumberFrames(ext, path = './'):
  fnames = glob('*.' + ext)
  return len(set(map(GetFrameId, fnames)))

def CombineNetcdfTiles(case, path = './'):
  nblocks = GetNumberBlocks('nc', path)

def ReadParticleAscii(case, name, path = './'):
  fnames = glob(path + '%s.%s.block*.?????.pat' % (case, name))
  fnames = sorted(fnames, key = lambda x: (GetFrameId(x), GetBlockId(x)))

  nblocks = GetNumberBlocks('pat', path)
  nframes = GetNumberFrames('pat', path)

  data, time = [], []

  for fname in fnames:
    fid = GetFrameId(fname)
    data.append(genfromtxt(fname))
    with open(fname, 'r') as ff:
      header = ff.readline()
    timestr = re.findall('time=\d+\.\d+e[+-]\d+', header)[0]
    time.append(float(timestr[5:]))

  for i in range(nframes):
    data[i*nblocks] = vstack(data[i*nblocks:(i+1)*nblocks])

  return data[::nblocks], time[::nblocks]
