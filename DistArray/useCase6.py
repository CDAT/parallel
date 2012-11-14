##UVCDAT use case 6
##Compute hybrid to std plevs


#!/usr/bin/env python

"""
From Michael Wehner, modified to run in parallel by Alex Pletzer (Tech-X)

To run this script, cpopy the following data:
hopper.nersc.gic:/project/projectdirs/vacet/Wehner/uvcdat_tests/

PS_cam5_1_amip_run2.cam2.h0.1994.nc
V_cam5_1_amip_run2.cam2.h0.1994.nc
cam5_1_amip_run2.cam2.h0.2002-09.nc

To execute the script:

time mpiexec -n 12 useCase6.py ./ V_cam5_1_amip_run2.cam2.h0.1994.nc V cam5_1_amip_run2.cam2.h0.2002-09.nc PS_cam5_1_amip_run2.cam2.h0.1994.nc 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10

On multipole.txcorp.com 
Linux multipole.txcorp.com 2.6.32-220.17.1.el6.x86_64 #1 SMP Wed May 16 00:01:37 BST 2012 x86_64 x86_64 x86_64 GNU/Linux
has 8 procs AMD Opteron(tm) Processor 4130 2.6Ghz

1 proc: 
real	16m50.452s
user	15m49.897s
sys	0m38.921s

4 procs:
real	5m25.258s
user	20m0.773s
sys	0m49.133s

6 procs:
real	4m35.794s
user	22m10.534s
sys	3m44.607s

12 procs:
real	5m36.064s
user	24m55.288s
sys	9m30.538s

Note: there are 12 time steps. Parallelization is over time. There is a critical
region where the data are written to the Netcdf file in order. Further improvement
in scaling can be obtained by having each processor write to a different file.

"""


# Adapted for numpy/ma/cdms2 by convertcdms.py
import numpy.oldnumeric as Numeric
# Input arguments to the script

from mpi4py import MPI

import sys,numpy.ma as MA,string,cdms2 as cdms,cdtime,cdutil,numpy.oldnumeric as Numeric
from cdtime import reltime
# reset DefaultCalendar below
cdtime.DefaultCalendar=cdtime.NoLeapCalendar

cdms.setNetcdfShuffleFlag(0)
cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)

pth="."
infile="sample_case_6.nc"
varin="cl"
ABfile=infile
outfile='y'+infile
varout=varin
psfile=infile
levels=[1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10]

Avar='a'
Bvar='b'
P0var='variable_13'

varps='ps'
levels=MA.asarray(list(levels))*100.
# Define some neat function

def getMetadata(s,**args):
  return s.getAxisList(**args),s.attributes
##   return apply(s.getAxisList,(),**args),s.attributes


def putMetadata(meta,tv,**args):
  import cdms2 as cdms
  if cdms.isVariable(meta) : meta=getMetadata(meta,**args)
  return cdms.createVariable(tv,attributes=meta[1])


def asMA(b,**keyargs):
  '''Converts a list of MAs to an MA
  Version 1.0 - February 28th, 2001
  Author Paul Dubois (dubois1@llnl.gov)
  '''
  import numpy.ma as MA
  m = [MA.getmaskarray(x) for x in b]
  b = [MA.filled(x) for x in b]
  return MA.array(b, mask=m,**keyargs)


    

# Gets the coef from the AB file
#print 'Getting coeff'
f=cdms.open(ABfile)
A=f(Avar)
B=f(Bvar)
Po=f(P0var)
#print 'Po:',Po
# open and get the surface pressure
#print 'Getting surface pressure'
fps=cdms.open(psfile)

# get var on sigma var (for ex T)
#print 'Getting T !'
fi=cdms.open(infile)
tt=fi.dimensionobject('time')
if tt.calendar=='360_day':cdtime.DefaultCalendar=cdtime.Calendar360
if tt.calendar=='gregorian':cdtime.DefaultCalendar=cdtime.MixedCalendar
if tt.calendar=='365_day':cdtime.DefaultCalendar=cdtime.NoLeapCalendar
if tt.calendar=='noleap':cdtime.DefaultCalendar=cdtime.NoLeapCalendar
if tt.calendar=='proleptic_gregorian':cdtime.DefaultCalendar=cdtime.GregorianCalendar
if tt.calendar=='standard':cdtime.DefaultCalendar=cdtime.StandardCalendar

v=fi.variables[varin]
meta=getMetadata(v)
tim=v.getTime()
u=tim.units

nprocs = MPI.COMM_WORLD.Get_size()
myPe = MPI.COMM_WORLD.Get_rank()
ranks = [rk%nprocs for rk in range(len(tim))]
outs = [None for rk in range(len(tim))]
if myPe == 0:
  print 'ranks = ', ranks

itim=0
import MV2

axes=v.getAxisList()
opened = False
for itim,tm in enumerate(tim): # loop over time
 
    if ranks[itim] == myPe:

      T=v(time=slice(itim,itim+1))
      Ps=fps('ps',time=slice(itim,itim+1))

      # create the pressure field
      P=cdutil.reconstructPressureFromHybrid(Ps,A,B,Po)

      # interpolate
      print '[%d] Interpolating at time index %d' % (myPe, itim)
      out=cdutil.logLinearInterpolation(T,P,levels)
      out.info()
      ## print '[%d] Done!' % myPe
      ## sh=list(out.shape)
      ## sh.insert(0,1)
      ## out=MA.reshape(out,tuple(sh))
      ## t=tim.subAxis(itim,itim+1)
      ## xx=reltime(tim[itim],tim.units)
      ## t_new=xx.torel('days since 1800').value
      ## t[0]=t_new
      ## t.units='days since 1800'
      ## meta[0][0]=t
      ## levelsax=cdms.createAxis(levels/100.)
      ## levelsax.id='plev'
      ## levelsax.units='hPa'
      ## levelsax.designateLevel()
      ## meta[0][1]=levelsax
      ## out=putMetadata(meta,out)
      ## out.info()
      ## out.id=varout
      outs[itim]=out
for itim in range(len(tim)):
    out = outs[itim]
    if out is not None: 
        if itim==0:
            f=cdms.open(outfile,'w')
        else:
            f=cdms.open(outfile,'r+')
        print '[%d] Writing data for itim = %d...' % (myPe, itim)
        f.write(out,id=varin)
        print '[%d] Done!' % myPe
        f.close()
    
    # make sure the data are written in order
    MPI.COMM_WORLD.Barrier()
