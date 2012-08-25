from mpi4py import MPI
import numpy
comm = MPI.COMM_WORLD
import cdms2,MV2

def addOne(starts,lens,j):
    starts[j]+=1
    if starts[j]==lens[j]:
        starts[j]=0
        starts=addOne(starts,lens,j+1)
    return starts

class Decompose(object):
    def __init__(self,V,axes=0):
        self.axes=[]
        for a in V.getAxisList():
            self.axes.append(a.clone())
        self.Aproc = comm.Get_rank() # This processor id
        self.Nproc = comm.Get_size() # N processors total
        print "OK:",self.Aproc,self.Nproc
        
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
        strides = int(numpy.power(self.Nproc,1./len(axes)))
        ls = []
        N=self.Nproc
        for i in range(len(axes)-1):
            ls.append(strides)
            N/=strides
            #print N
        ls.append(N)
        N=1
        for n in ls: N*=n
        if N!=self.Nproc:
            ls[-1]-=1
            #raise Exception,"Can't decompose %i proc over %i dims" % (self.Nproc,len(axes))

        kw={}
        starts=[0,]*len(axes)
        i=0
        while i<self.Aproc:
            starts=addOne(starts,ls,0)
            i+=1


        #print 'Starts:',Aproc,starts
        j=0
        kw={}
        self.ids = []
        for i in range(len(V.shape)):
            if i in Is:
                self.ids.append(V.getAxis(i).id)
                kw[self.ids[-1]]=slice(starts[j],None,ls[j])
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
            
        if self.Aproc == root:
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
            for i in range(self.Nproc):
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
                
