#! /usr/bin/env python2.7
from glob import glob
from numpy import genfromtxt, zeros, vstack 
from subprocess import check_call
import re, os, getopt, sys, operator, pickle

def List2str(s, format = None):
  if format == None:
    return reduce(operator.add, [' ' + str(x) for x in s[1:]], str(s[0]))
  else:
    return reduce(operator.add, [format % x for x in s[1:]], format % s[0])

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
  fnames = glob(path + '/*.[0-9][0-9][0-9][0-9][0-9].' + ext)
  return len(set(map(GetBlockId, fnames)))

def GetNumberFrames(ext, path = './'):
  fnames = glob(path + '/*.?????.' + ext)
  return len(set(map(GetFrameId, fnames)))

def GetNcVariable(fname, vname, path = './'):
  data = Dataset(fname, 'r')
  return squeeze(data[vname][:])

def CombineNcFiles(case, path = './', out = ''):
  print 'combining netcdf files ...'
  if out == '':
    outfile = case
  else:
    outfile = case + '-' + out
  nblocks = GetNumberBlocks('nc', path)
  for i in range(nblocks):
    print 'processing block %d ...' % i
    files = path + '/' + '%s.block%d.*.*.nc' % (case, i)
    target = outfile + '.nc.%04d' % i
    check_call('ncrcat -h %s -o %s' % (files, target), shell = True)
    check_call('ncatted -O -a %s,%s,%c,%c,%d %s' % ('NumFilesInSet', 'global', 'c', 'i',
              nblocks, target), shell = True)

  #if not os.path.exists('mppnccombine'):
  #  subprocess.call('gcc -O -o mppnccombine mppnccombine.c -lnetcdf', shell = True)

  check_call('./mppnccombine %s.nc' % outfile, shell = True)

  for f in glob(outfile + '.nc.????'):
    os.remove(f)
  #for f in glob(path + '/' + case + '.block*.nc'):
  #  os.remove(f)

  print 'Output file writted in ', outfile + '.nc'

def ReadParticleAscii(case, name, path = './', fs = None, fe = None):
  fnames = glob(path + '/' + '%s.%s.block*.?????.pat' % (case, name))
  fnames = sorted(fnames, key = lambda x: (GetFrameId(x), GetBlockId(x)))

  nblocks = GetNumberBlocks('pat', path)
  nframes = GetNumberFrames('pat', path)
  print nblocks, nframes
  if fs == None: fs = 0
  if fe == None: fe = nframes

  data, time = [], []

  for i, fname in enumerate(fnames[fs*nblocks:fe*nblocks]):
    if i%nblocks == 0:
      print 'processing frame %d ...' % (fs+i/nblocks,)
    fid = GetFrameId(fname)
    data.append(genfromtxt(fname))
    with open(fname, 'r') as ff:
      header = ff.readline()
    timestr = re.findall('time=\d+\.\d+e[+-]\d+', header)[0]
    time.append(float(timestr[5:]))

  if nblocks > 1:
    for i in range(fe - fs):
      data[i*nblocks] = vstack(data[i*nblocks:(i+1)*nblocks])

  return data[::nblocks], time[::nblocks]

def CombinePatFiles(case, name, path = './', out = ''):
  print 'combining pat files ...'
  if out == '':
    outfile = case
  else:
    outfile = case + '-' + out
  data, time = ReadParticleAscii(case, name, path)
  with open(outfile + '.pat.pk', 'wb') as file:
    pickle.dump(data, file)
    pickle.dump(time, file)

def ReadCombinedPatFiles(patfile):
  with open(patfile, 'rb') as file:
    data = pickle.load(file)
    time = pickle.load(file)
  return data, time

if __name__ == '__main__':
  opts, args = getopt.getopt(sys.argv[1:], 'c:d:no:p:', ['case', 'dir', 'netcdf',
    'output', 'particle'])

  task_nc_combine   = False
  task_pat_combine  = False

  case_name = ''
  path_name = './'
  out_name  = ''
  particle_name = ''

  for opt, arg in opts:
    if opt in ['-c', '--case']:
      case_name = arg
    elif opt in ['-d', '--dir']:
      path_name = arg
      out_name = arg
    elif opt in ['-n', '--netcdf']:
      task_nc_combine = True
    elif opt in ['-o', '--output']:
      out_name = arg
    elif opt in ['-p', '--pat']:
      task_pat_combine = True
      particle_name = arg
    else:
      assert False, 'unhandled option'

  if task_nc_combine:
    CombineNcFiles(case_name, path = path_name, out = out_name)

  if task_pat_combine:
    CombinePatFiles(case_name, particle_name, path = path_name, out = out_name)
