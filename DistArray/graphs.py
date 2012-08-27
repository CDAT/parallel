import sys
import numpy
import cdms2
import MV2

test = sys.argv[1]

cpus = [1,2,3,4,6,8]
methods = ["stagger","blocks"]
speed=numpy.zeros((10,len(methods),len(cpus),3))
machines=[]
os=[]
ver=[]
for i,n in enumerate(cpus):
    for j,method in enumerate(methods):
        files = os.listdir(".")
        for fnm in files:
            if not (fnm[:7]=="timing_" and fnm.find("_%s_" % test)>-1):
                continue
            sp = fnm.split("_")
            
            f=open(fnm)
            k=0
            while True:
                try:
                    rd = float(f.readline().split()[1])
                    ex = float(f.readline().split()[1])
                    speed[j,i,0] += rd
                    speed[j,i,1] += ex - rd
                    speed[j,i,2] += ex 
                    k+=1
                except:
                    break
            speed[j,i,0]/=k
            speed[j,i,1]/=k
            speed[j,i,2]/=k
cpu = cdms2.createAxis(range(len(cpus)))
cpu = cdms2.createAxis(cpus)
cpu.id = "Ncpus"
mthd =cdms2.createAxis(range(len(methods)))
mthd.id='Method'


import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=3,ncols=speed.shape[0])

print ax.shape


ax[0,0].set_ylabel('Read Time (s)')
ax[1,0].set_ylabel('Compute Time (s)')
ax[2,0].set_ylabel('Total Time (s)')
for i in range(speed.shape[0]):
    ax[0,i].set_title(methods[i])
    for j in range(3):
        ax[j,i].plot(cpus,speed[i,:,j],'-o',ms=5,mfc='orange')
    #ax[0,i].set_xlabel('NCPUS')
    ax[-1,i].set_xlabel('NCPUS')
plt.savefig("DJF.png")
plt.show()

