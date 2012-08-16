import cdms2
import cdutil
import time
import vertical
import numpy

f=cdms2.open("hus_3hrMlev_HadGEM2-A_tamip200907_r1i1p1_200907150130-200907192230.nc")

b=f("b")
orog = f("orog")
lev = b.getLevel()+(b-b)

start=time.time()
z = cdutil.vertical.reconstructPressureFromHybrid(orog,lev,b,1.)
duration = time.time()-start
print z.shape,duration*1000,"ms"

sh = list(orog.shape)
sh.insert(1,len(b))
dummy= numpy.zeros(sh)
orog=orog.astype("d")
print orog.dtype,lev.dtype,b.dtype,dummy.dtype
start=time.time()
z2 = vertical.numbaReconstruct4D(orog.filled(),lev.filled(),b.filled(),1.,dummy)
#z2 = vertical.reconstructPressureFromHybrid(orog,lev,b,1.)
duration = time.time()-start
print z2.shape,duration*1000,"ms"
