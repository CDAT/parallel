# Adapted for numpy/ma/cdms2 by convertcdms.py
import MV2
import numpy
import numba
from numba.decorators import jit as jit

def recontruct4DCore(ps,A,B,P0,dummy):
    n1,n2,n3,n4 = dummy.shape
    p = numpy.zeros_like(dummy)
    for j in range(n2):
        Afact=A[j]*P0
        for i in range(n1):
            for k in range(n3):
                for l in range(n4):
                    p[i,j,k,l]=ps[i,k,l]*B[j]+Afact
    return p

numbaReconstruct4D = jit(ret_type = numba.d[:,:,:,:], arg_types=[numba.d[:,:,:],numba.d[:],numba.d[:],numba.d,numba.d[:,:,:,:]])(recontruct4DCore)

                                          
def reconstructPressureFromHybrid(ps,A,B,Po):
    """
    Reconstruct the Pressure field on sigma levels, from the surface pressure
    
    Input
    Ps   : Surface pressure
    A,B,Po: Hybrid Convertion Coefficients, such as: p=B.ps+A.Po
    Ps: surface pressure
    B,A are 1D : sigma levels
    Po and Ps must have same units
    
    Output
    Pressure field
    Such as P=B*Ps+A*Po

    Example
    P=reconstructPressureFromHybrid(ps,A,B,Po)
    """
    # Compute the pressure for the sigma levels
    ax1 = ps.getAxisList()
    ax2 = A.getAxisList()
    sh = list(ps.shape)
    sh.insert(1,len(A))
    dummy= numpy.zeros(sh)


    p = MV2.array(numbaReconstruct4D(ps.filled(),A.filled(),B.filled(),Po,dummy))
    p.setAxisList((ax1[0],ax2,ax1[1],ax1[3]))

    p.id='P'
    try:
      p.units=ps.units
    except:
      pass
    t=ps.getTime()
    if not t is None:
      p=p(order='tz...')
    else:
     p=p(order='z...')
    return p
    

## def LinearCore4d(A,P,levels):
##     nlev=levels.shape[0]
##     nsigma,nt,ny,nx=A.shape
##     out = numpy.zeros((nlev,nt,ny,nx))
##     for ilev in range(nlev): # loop through pressure levels
##         lev=levels[ilev] # get value for the level
##         Pabv=-numpy.ones((nt,ny,nx))
##         Aabv=-numpy.ones((nt,ny,nx)) # Array on sigma level Above
##         Abel=-numpy.ones((nt,ny,nx)) # Array on sigma level Below
##         Pbel=-numpy.ones((nt,ny,nx)) # Pressure on sigma level Below
##         Pabv=-numpy.ones((nt,ny,nx)) # Pressure on sigma level Above
##         Peq=-1.e20*numpy.ones((nt,ny,nx)) # Area where Pressure == levels
##         for i in range(1,nsigma): # loop from second sigma level to last one
##             for j in range(nt):
##                 for k in range(ny):
##                     for l in range(nx):
##                         if P[i,j,k,l]>=lev and P[i-1,j,k,l]<=lev:
##                             Pabv[j,k,l] = P[i,j,k,l]
##                             Aabv[j,k,l] = A[i,j,k,l]
##                             Pbel[j,k,l] = P[i-1,j,k,l]
##                             Abel[j,k,l] = A[i-1,j,k,l]
##                         if P[i,j,k,l] == lev:
##                             Peq[j,k,l] = A[i,j,k,l]

##         for j in range(nt):
##             for k in range(ny):
##                 for l in range(nx):
##                     if Pbel == -1.:
##                         out[ilev,j,k,l] = 1.e20
##                     else:
##                         if Peq == -1.e20:
##                             out[ilev,j,k,l] = (lev-Pbel[j,k,l])/(Pabv[j,k,l]-Pbel[j,k,l])*(Aabv[j,k,l]-Abel[j,k,l])+Abel[j,k,l]
##                         else:
##                             out[ilev,j,k,l] = Peq[j,k,l]
                    
##         return out
## numbaLinearCore4D = numba.decorators.jit(ret_type = numba.d[:,:,:,:], arg_types=[numba.d[:,:,:],numba.d[:,:,:],numba.d[:]])(LinearCore4d)
    
## def LinearInterpolation(A,P,levels):
##     """
##     linear interpolation
##     to convert a field from sigma levels to pressure levels
##     Value below surface are masked
    
##     Input
##     A :    array on sigma levels
##     P :    pressure field from TOP (level 0) to BOTTOM (last level)
##     levels : pressure levels to interplate to (same units as P), default levels are:[100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000]

##     P and levels must have same units

##     Output
##     array on pressure levels (levels)
    
##     Examples:
##     A=logLinearInterpolation(A,P),levels=[100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000])
##     """
    
##     try:
##         nlev=len(levels)  # Number of pressure levels
##     except:
##         nlev=1  # if only one level len(levels) would breaks
##         levels=[levels,]
##     order=A.getOrder()
##     A=A(order='z...')
##     P=P(order='z...')
##     sh=list(P.shape)
##     nsigma=sh[0] #number of sigma levels
##     sh[0]=nlev
##     t=MV2.zeros(sh,typecode=MV2.float32)
##     sh2=P[0].shape
        
##     ax=A.getAxisList()
##     autobnds=cdms2.getAutoBounds()
##     cdms2.setAutoBounds('off')
##     lvl=cdms2.createAxis(MV2.array(levels).filled())
##     cdms2.setAutoBounds(autobnds)
##     try:
##         lvl.units=P.units
##     except:
##         pass
##     lvl.id='plev'
    
##     try:
##       t.units=P.units
##     except:
##       pass
  
##     ax[0]=lvl
##     t.setAxisList(ax)
##     t.id=A.id
##     for att in A.listattributes():
##         setattr(t,att,getattr(A,att))
##     return t(order=order)
    
## sigma2Pressure=logLinearInterpolation
