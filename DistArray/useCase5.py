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

D = Decompose(V,axes=['x','y'])

sub = V(D())
print "Decomped",D.Aproc,sub.shape,sub.count()
cdutil.times.setTimeBoundsMonthly(sub)
sub=cdutil.DJF.departures(sub)
print "DJFs",D.Aproc,sub.shape,sub.count()

print "Gathering"

out = D.gather(sub)

if D.Aproc == 0 :
    print out.shape

if D.Aproc == 0:
    import vcs
    x=vcs.init()
    x.plot(out)
    raw_input()
#print D.Aproc,out
