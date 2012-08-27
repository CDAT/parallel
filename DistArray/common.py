from mpi4py import MPI
import numpy
comm = MPI.COMM_WORLD
import cdms2,MV2
import time
import sys
import os



Aproc = comm.Get_rank() # This processor id
Nproc = comm.Get_size() # N processors total

def addOne(starts,lens,j):
    starts[j]+=1
    if starts[j]==lens[j]:
        starts[j]=0
        starts=addOne(starts,lens,j+1)
    return starts

class Decompose(object):
    def __init__(self,V,axes=0,method='stagger'):
        self.axes=[]
        for a in V.getAxisList():
            self.axes.append(a.clone())
        self.n = comm.Get_rank() # This processor id
        self.N = comm.Get_size() # N processors total
        
        if isinstance(axes,int):
            axes=[axes,]

        Ns = []
        Is = []
        for a in axes:
            if isinstance(a,int):
                Ns.append(V.shape[a])
                Is.append(a)
            elif isinstance(a,str):
                if not a in ['x','y','z','t']:
                    for i,nm in enumerate(V.getAxisIds()):
                        if a==nm:
                            Ns.append(V.shape[i])
                            Is.append(i)
                            break
                else:
                    for i,ax in enumerate(V.getAxisList()):
                        if ((a=='t' and ax.isTime()) or
                            (a=='z' and ax.isLevel()) or
                            (a=='y' and ax.isLatitude()) or
                            (a=='x' and ax.isLongitude())
                            ):
                            Ns.append(len(ax))
                            Is.append(i)

        Ns.sort(key=lambda Is:Is)
        Is.sort()
        print "Ns,Is:",Ns,Is
        strides = int(numpy.power(self.N,1./len(axes)))
        ls = []
        N=self.N
        for i in range(len(axes)-1):
            ls.append(strides)
            N/=strides
            #print N
        ls.append(N)
        N=1
        for n in ls: N*=n
        if N!=self.N:
            ls[-1]-=1
            #raise Exception,"Can't decompose %i proc over %i dims" % (self.N,len(axes))

        kw={}
        starts=[0,]*len(axes)
        i=0
        while i<self.n:
            starts=addOne(starts,ls,0)
            i+=1


        #print 'Starts:',n,starts
        j=0
        kw={}
        self.ids = []
        for i in range(len(V.shape)):
            if i in Is:
                ax = V.getAxis(i)
                self.ids.append(ax.id)
                if method=='stagger':
                    kw[self.ids[-1]]=slice(starts[j],None,ls[j])
                elif method=='blocks':
                    l = len(ax)/ls[j]
                    kw[self.ids[-1]]=slice(starts[j]*l,(starts[j]+1)*l)
                else:
                    raise Exception, "Unknown method: %s" % method
                j+=1

        self.selector = kw
        self.pos = Is
        return

    def __call__(self,slab=None):
        S=cdms2.selectors.Selector(**self.selector)
        if slab is None:
            return S
        else:
            return slab(S)

    def getSlices(self):
        slices = []
        for i in self.ids:
            slices.append(self.selector[i])
        return slices
        
    def gather(self,array,pos=None,root=0):

        lst = comm.gather(array)
        lst2 = comm.gather(self)
        if pos is None:
            pos=self.pos
            
        if self.n == root:
            sh=[]
            axes=[]
            j=0
            for i,a in enumerate(self.axes):
                if i in pos:
                    sh.append(len(a))
                    axes.append(self.axes[i])
                else:
                    sh.append(array.shape[i])
                    axes.append(array.getAxis(i))
                
                    
            out=numpy.zeros(sh)
            for i in range(self.N):
                slices=lst2[i].getSlices()
                print i,slices
                Slices=[]
                k=0
                for j in range(len(self.axes)):
                    print k,j,pos
                    if j in pos:
                        Slices.append(slices[k])
                        k+=1
                    else:
                        Slices.append(slice(0,None))
                out[Slices]=lst[i]
            out=MV2.array(out)
            out.setAxisList(axes)
            out.id=array.id
            return out
                
start=time.time()

def writeTime(method,txt):
    end=time.time()
    args = (sys.argv[0][:-3],)+os.uname()[:3]+(Nproc,method)
    print args
    fnm="timing_%s_%s_%s_%s_%i_%s.txt" % args
    if os.path.exists(fnm):
        f=open(fnm,"r+")
        f.read()
    else:
        f=open(fnm,"w")
    print >>f, txt,end-start
    print txt,end-start
    f.close()
