##UVCDAT Use case 5
## DJF departures

## We will split the data over lat/lon

from common import *

import cdutil
import cdms2
import sys
fnm,vr = sys.prefix+'/sample_data/clt.nc','clt'
fnm,vr = 'tls_both_ccsm4.hist_rcp85.r1i1p1.mo.nc','eqmsu_tls'
f=cdms2.open(fnm)

V=f[vr]

if len(sys.argv)>1:
    method=sys.argv[1]
else:
    method='blocks'
    
D = Decompose(V,axes=['x','y'],method=method)

sub = V(D())
if D.n == 0 :
    writeTime(method,"read:")
print "Decomped",D.n,sub.shape,sub.count(),sub.getLatitude()[:5],sub.getLongitude()[:5]
#print "Sel:",D()
cdutil.times.setTimeBoundsMonthly(sub)
sub=cdutil.DJF.departures(sub)
if D.n == 0 :
    writeTime(method,"exec:")
print "DJFs",D.n,sub.shape,sub.count()

print "Gathering"

out = D.gather(sub)

if D.n == 0 :
    print out.shape
    writeTime(method,"gather:")

## if D.n == 0:
##     import vcs
##     x=vcs.init()
##     x.plot(out)
##     raw_input()
## #print D.n,out
