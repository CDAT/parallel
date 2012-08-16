import cdms2
import MV2
import distarray
import numpy
import cdutil
from mpi4py import MPI

# open file and read data
clt = cdms2.open('/usr/local/uvcdat/1.0.0/sample_data/clt.nc', 'r')['clt']

def getDomainDecomp(nprocs,sizes):
    n2=nprocs//2
    res0 = [1,nprocs]
    res = res0[:]
    k=2
    while k <= n2:
        if (res0[0]*k) * (res0[1]//k) == nprocs:
            res[0]=res0[0]*k
            res[1] = res0[1]//k
        k += 1
    return res

nLat,nLon = clt.shape[1:]
rk = MPI.COMM_WORLD.Get_rank()
sz = MPI.COMM_WORLD.Get_size()

npLat,npLon = getDomainDecomp(sz,(nLat,nLon))

if rk == 0:
    value=0
    cdms2.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
    cdms2.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
    cdms2.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included
    print "nlon,nlat:",nLat,nLon
print rk,"nps:",npLat,npLon

nlat,nlon = nLat/npLat , nLon/npLon
ipLat = rk // npLon
ipLon = rk % npLon

iLatBeg , iLatEnd = ipLat * nlat, (ipLat+1) * nlat
iLonBeg , iLonEnd = ipLon * nlon, (ipLon+1) * nlon

daclt = clt[:,iLatBeg:iLatEnd,iLonBeg:iLonEnd]
cdutil.setTimeBoundsMonthly(daclt)
print rk,daclt.shape
mp = cdutil.averager(daclt,axis='t')


if rk==0:
    print "Gathering results"
lst = MPI.COMM_WORLD.gather(mp,root=0)
if rk==0:
    print "Gathered"
    out = numpy.zeros(clt.shape[1:],clt.dtype)
    for proc in range(sz):
        print len(lst),proc,lst[proc].shape
    for proc in range(sz):
        ipLat = proc // npLon
        ipLon = proc % npLon
        iLatBeg , iLatEnd = ipLat * nlat, (ipLat+1) * nlat
        iLonBeg , iLonEnd = ipLon * nlon, (ipLon+1) * nlon
        out[iLatBeg:iLatEnd,iLonBeg:iLonEnd] = lst[proc]

    lat = clt.getLatitude()
    lon = clt.getLongitude()

    out= MV2.array(out)
    print out.shape
    out.setAxis(0,lat)
    out.setAxis(1,lon)
    f=cdms2.open("map.nc","w")
    f.write(out,id=clt.id)
    f.close()
    
    import vcs
    x=vcs.init()
    x.plot(out)
    raw_input()
    
