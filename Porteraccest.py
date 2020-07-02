import sys

def function(x,a,b):
    return a * x + b

def porteracctest(popt,predfile,confstartloc,q3startloc):
    predf = open(predfile,'r').readlines()
    currline = 0
    diffs = []
    estaccs = []
    predseqs = []
    q3s = []
    for i in range(1,len(truef),2):
        prevline = currline
        currline += 1
        while 'Query_name' not in predf[currline] and currline < len(predf) - 1 and 'Query served' not in predf[currline]:
            currline += 1
        conf = getconf(predf[prevline:currline],confstartloc)
        a = 0.06717725 
        b = 0.3950873643451706
        estacc = function(conf,a,b)
        estaccs.append(estacc)
    estpred = [(estaccs[i],predseqs[i],q3s[i]) for i in range(len(estaccs))]
    #estpred = [(estaccs[i],predseqs[i]) for i in range(len(estaccs))]
    return estpred

