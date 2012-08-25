from mpi4py import MPI
    
def sliding(foo,chunkSize,*args):
    """ Applies foo over nProc processors with a sliding chunk of size "chunkSize", over axis "chunkAxis"
    """

    N = 0
    for a in args:
        if isinstance(a,numpy.ndarray):
            N=a.shape[0]
            

    Aproc = MPI.COMM_WORLD.Get_rank() # This processor id
    Nproc = MPI.COMM_WORLD.Get_size() # N processor total

    ## Ok now 
    
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
