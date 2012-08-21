import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import cdms2 as cdms
#from RMT import get_grid
import cdutil
from mpi4py import MPI
import os

def get_data():
    f=cdms.open("/export/marvel1/RandomMatrix/TLS/ccsm4/tls_both_ccsm4.hist_rcp85.r1i1p1.mo.nc")

    da=cdutil.ANNUALCYCLE.departures(f("eqmsu_tls"))
    GRID=cdms.createUniformGrid(-90,9,20,-180,18,20)

    return da.regrid(GRID,regridTool='esmf')

def getgrid(grid="coarse"):
    file1 ='/usr/local/uvcdat/2012-06-11/sample_data/hadcrut2_sample.nc'
    a = cdms.open(file1) 
    data=a('temanom')
    if grid=="coarse":
        
        LATS=data.getLatitude()[::4]
        ub=data.getLatitude().getBounds()[1::4][:,0]
        lb=data.getLatitude().getBounds()[::4][:,1]
        latbounds=np.array(zip(ub,lb))
        
        LONS=data.getLongitude()[::4]
        ub=data.getLongitude().getBounds()[1::4][:,0]
        lb=data.getLongitude().getBounds()[::4][:,1]
        lonbounds=np.array(zip(ub,lb))
    else:
        
        LATS=data.getLatitude()[::2]
        ub=data.getLatitude().getBounds()[1::2][:,0]
        lb=data.getLatitude().getBounds()[::2][:,1]
        latbounds=np.array(zip(ub,lb))
        
        LONS=data.getLongitude()[::2]
        ub=data.getLongitude().getBounds()[1::2][:,0]
        lb=data.getLongitude().getBounds()[::2][:,1]
        lonbounds=np.array(zip(ub,lb))
    return cdms.createGenericGrid(LATS,LONS,latBounds=latbounds,lonBounds=lonbounds)



def IPR(vec):
    return np.sum(vec**4)

def PR(vec):
    return 1./(len(vec)*IPR(vec))

def pr_spectrum(vals,vecs):
    vals=np.real(vals)
    ncomp = float(len(vals))
    I=np.argsort(vals)[::-1]
    A=vecs[:,I]#sort in descending order of eigenvalue
    return 1./(ncomp*np.sum(A**4,axis=0))
    
def chunky(test,chunk_size=240):
   
    nt=test.shape[0]
    rk = MPI.COMM_WORLD.Get_rank()
    sz=MPI.COMM_WORLD.Get_size()
    
    ns=test.shape[1]*test.shape[2]
    #particip=np.zeros((nt,ns))
    particip=[]
    #eigenvalues=np.zeros((nt,ns))
    eigenvalues=[]
    i=0
    
    while i< nt-chunk_size-rk:
       
        start=i+rk
        stop=start+chunk_size
        #print "start is "+str(start)
        #print "stop is "+str(stop)
        A=test[start:stop]
        A=A.reshape(A.shape[0],A.shape[1]*A.shape[2])
        C=np.corrcoef(A.T)
        if C.shape != (162,162):
            print "rk is "+str(rk)
            print"++++++++++++++++"
            print "A.T.shape is ",A.T.shape
            print "i is:",i,"start",start,"end",stop
            print "********************"
            
            print "C.shape is "+str(C.shape)
        vals,vecs=la.eig(C)
        
        #particip[start]=pr_spectrum(vals,vecs)
        #eigenvalues[start]=np.real(vals)
        eigenvalues.append(np.real(vals))
        particip.append(pr_spectrum(vals,vecs))
            
        i+=sz
    if rk==0:
        print "Done, converting to numpy and going back"
    eigenvalues=np.array(eigenvalues)
    particip=np.array(particip)
    return eigenvalues,particip
     

chunk_size=240
test=f=cdms.open("/export/marvel1/RandomMatrix/TLS/ccsm4/tls_both_ccsm4.hist_rcp85.r1i1p1.mo.REGRID.nc")['TLS']
nt,n1,n2=test.shape
ns=n1*n2
#EV,D=chunky(test)
eigenvalues,particip=chunky(test,chunk_size=chunk_size)
rk = MPI.COMM_WORLD.Get_rank()
sz=MPI.COMM_WORLD.Get_size()
if rk==0:
    print "Gathering now"
print "At the end:",rk,eigenvalues.shape

ge = MPI.COMM_WORLD.gather(eigenvalues,root=0)
gp = MPI.COMM_WORLD.gather(particip,root=0)

if rk==0:
    print "gathered:",len(ge),ge[0].shape

if rk==0:
    E = np.zeros((nt-chunk_size,ns))
    P = np.zeros((nt-chunk_size,ns))
    i=0
    while i < ge[0].shape[0]:
        start=i
        for j in range(sz):
            E[i*sz+j]=ge[j][i]
            P[i*sz+j]=gp[j][i]
        i+=1
        


    #eigenvalues[start]=MPI.COMM_WORLD.gather(eigenvalues,root=0)
    #particip[start]=MPI.COMM_WORLD.gather(particip,root=0)


    import pickle
    f=open("eigenvectors_particip.dat","w")
    pickle.dump([E,P],f)
    f.close()
