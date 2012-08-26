import sys,os

test = sys.argv[1]

for n in [1,2,3,4,6,8]:
    for method in ["stagger","blocks"]:
        cmd = "mpiexec -n %i python %s %s" % (n,test,method)
        for i in range(5):
            print "%i, cmd: %s" % (i,cmd)
            ln = os.popen(cmd).readlines()
