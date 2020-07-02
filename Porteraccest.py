import sys

def function(x,a,b):
    return a * x + b

def getconf(lines,startloc,returnarray=False):
    conf = 0
    count = 0
    outintline = []
    for i in range(startloc,len(lines),7):
        line = lines[i].strip()
        intline = [int(l) for l in line.strip()]
        outintline += intline
        conf += sum(intline)
        count += len(line)
    if returnarray == False:
        return conf / count
    else:
        return conf/count,outintline

def porteracctest(predfile,confstartloc,q3startloc,a,b):
    predf = open(predfile,'r').readlines()
    currline = 0
    diffs = []
    estaccs = []
    predseqs = []
    q3s = []
    prevline = currline
    currline += 1
    while 'Query_name' not in predf[currline] and currline < len(predf) - 1 and 'Query served' not in predf[currline]:
        currline += 1
    conf = getconf(predf[prevline:currline],confstartloc)
    estacc = function(conf,a,b)
    return estacc

predfile = sys.argv[1]
if sys.argv[2] == '8state':
    confstartloc = 7
    q3startloc = 6
    a = 0.06717725 
    b = 0.3950873643451706
elif sys.argv[2] == '3state':
    confstartloc = 5
    q3startloc = 4
    a = 0.05381087 
    b = 0.5139789336088341

print('estimated accuracy: {:.03f}'.format(porteracctest(predfile,confstartloc,q3startloc,a,b)))
