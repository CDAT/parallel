import sys
import numpy
import cdms2
import MV2
import os

test = sys.argv[1]

#cpus = [1,2,3,4,6,8]
#methods = ["stagger","blocks"]
machines=[]
Os=[]
ver=[]
methods=[]
cpus=[]
Files = os.listdir(".")
files=[]


def parseFileName(fnm):
    sp=fnm.split("_")
    T = sp[1]
    O = sp[2]
    M = sp[3].split(".")[0]
    V = sp[4]
    n=int(sp[-2])
    m=sp[-1].split(".")[0]
    return T,O,M,V,n,m

    
for fnm in Files:
    if not (fnm[:7]=="timing_" and fnm.find("_%s_" % test)>-1):
        continue
    T,O,M,V,n,m = parseFileName(fnm)
    print O,M,V
    if not M in machines:
        machines.append(M)
        Os.append(O)
        ver.append(V)
    if not n in cpus:
        cpus.append(n)
    if not m in methods:
        methods.append(m)
    files.append(fnm)

speed=numpy.zeros((len(machines),len(methods),len(cpus),4))
cpus.sort()

for fnm in files:
    T,O,M,V,n,m = parseFileName(fnm)
    i=machines.index(M)
    j=methods.index(m)
    k=cpus.index(n)
    f=open(fnm)
    N=0
    while True:
        try:
            rd = float(f.readline().split()[1])
            ex = float(f.readline().split()[1])
            ga = float(f.readline().split()[1])
            speed[i,j,k,0] += rd
            speed[i,j,k,1] += ex - rd
            speed[i,j,k,2] += ga - ex 
            speed[i,j,k,3] += ga 
            N+=1
        except:
            break
    speed[i,j,k,0]/=N
    speed[i,j,k,1]/=N
    speed[i,j,k,2]/=N
    speed[i,j,k,3]/=N

speed=numpy.ma.masked_equal(speed,0)
#speed[:,:,:,3]= speed[:,:,:,3] / speed[:,:,0,3]

cpu = cdms2.createAxis(cpus)
cpu.id = "Ncpus"
mthd =cdms2.createAxis(range(len(methods)))
mthd.id='Method'


import matplotlib.pyplot as plt


fig, ax = plt.subplots(nrows=4,ncols=speed.shape[1])
fig.set_size_inches(8.,11.5)
print ax.shape

plt.subplots_adjust(right=.78)

print cpus
print machines
print methods

desc = {"meryem":"Mac (Snw)\n16 CPU 16Gb RAM",
        "sofia":"Mac (Snw)\n8CPU 8Gb RAM",
        "crunchy":"Linux (RH6)\n16CPU 96Gb RAM*\n*:Users intense"
        }
ax[0,0].set_ylabel('Read Time (s)')
ax[1,0].set_ylabel('Compute Time (s)')
ax[2,0].set_ylabel('Gather Time (s)')
ax[3,0].set_ylabel('Total Time (s)')

lgd = []
for i in range(speed.shape[1]):
    ax[0,i].set_title(methods[i])
    for m in range(speed.shape[0]):
        for j in range(4):
            ax[2,i].set_ylim([0,2])
            kw = {'linewidth':1.5,
                  'ms':5,
                  'mfc':'orange'
                  }
            if j==0 and i==speed.shape[1]-1:
                kw['label']=machines[m]
            ax[j,i].plot(cpus,speed[m,i,:,j],'-o',**kw)
    #ax[0,i].set_xlabel('NCPUS')
    ax[-1,i].set_xlabel('NCPUS')
    
for m in machines:
    lgd.append(desc[m])

ax[0,-1].legend(lgd,shadow=True,fancybox=True,bbox_to_anchor=[1.72,1.15])#,font_size=1.)
L=ax[0,-1].get_legend()

ltext = L.get_texts()
ltext[0].set_fontsize(6.5)
F=L.get_frame()
F.set_facecolor("beige")
plt.savefig("DJF.png")
plt.show()

raw_input()
