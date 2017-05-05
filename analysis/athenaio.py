#! /usr/bin/env python2.7
from glob import glob
from numpy import genfromtxt, zeros, vstack 
import re, subprocess, os, getopt, sys

def GetBlockId(fname):
  for name in fname.split('.'):
    if name[:5] == 'block':
      return int(name[5:])
  assert False, 'Block id not found'

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
  assert False, 'Output id not found'

def GetNumberBlocks(ext, path = './'):
  fnames = glob('*.[0-9][0-9][0-9][0-9][0-9].' + ext)
  return len(set(map(GetBlockId, fnames)))

def GetNumberFrames(ext, path = './'):
  fnames = glob('?????.' + ext)
  return len(set(map(GetFrameId, fnames)))

def GetNcVariable(fname, vname, path = './'):
  data = Dataset(fname, 'r')
  return squeeze(data[vname][:])

def CombineNcFiles(case, path = './', out = ''):
  print 'combining netcdf files ...'
  nblocks = GetNumberBlocks('nc', path)
  for i in range(nblocks):
    print 'processing block %d ...' % i
    files = path + case + '.block%d.*.*.nc' % i
    target = case + '.nc.%04d' % i
    subprocess.call('ncrcat -h %s -o %s' % (files, target), shell = True)
    subprocess.call('ncatted -O -a %s,%s,%c,%c,%d %s' 
      % ('NumFilesInSet', 'global', 'c', 'i', nblocks, target),
      shell = True)

  #if not os.path.exists('mppnccombine'):
  #  subprocess.call('gcc -O -o mppnccombine mppnccombine.c -lnetcdf', shell = True)

  if out == '':
    subprocess.call('./mppnccombine %s.nc' % case, shell = True)
  else:
    subprocess.call('./mppnccombine %s-%s.nc' % (case, out), shell = True)

  for f in glob(case + '.nc.????'):
    os.remove(f)
  print 'finished.'

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

  if nblocks > 1:
    for i in range(nframes):
      data[i*nblocks] = vstack(data[i*nblocks:(i+1)*nblocks])

  return data[::nblocks], time[::nblocks]

if __name__ == '__main__':
  opts, args = getopt.getopt(sys.argv[1:], 'c:d:o', ['combine', 'dir', 'output'])

  task_nc_combine = False
  case_name = ''
  path_name = './'
  out_name = ''

  for opt, arg in opts:
    if opt in ('-c', '--combine'):
      task_nc_combine = True
      case_name = arg
    elif opt in ('-d', '--dir'):
      path_name = arg
    elif opt in ('-o', '--output'):
      out_name = arg
    else:
      assert False, 'unhandled option'

  if task_nc_combine:
    CombineNcFiles(case_name, path = path_name, out = out_name)
