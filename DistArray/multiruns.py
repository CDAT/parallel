import sys,os

test = sys.argv[1]

if len(sys.argv)>2:
    cpus=[]
    for a in sys.argv[2:]:
        cpus.append(int(a))
else:
    cpus = [1,2,3,4,6,8] 
for n in cpus:
    for method in ["stagger","blocks"]:
        cmd = "mpiexec -n %i python %s %s" % (n,test,method)
        for i in range(5):
            print "%i, cmd: %s" % (i,cmd)
            ln = os.popen(cmd).readlines()
