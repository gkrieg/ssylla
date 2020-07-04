import sys

def getestimate(fname):
    inf = open(fname).readlines()
    for line in inf:
        if 'JNETCONF' in line:
            confstr = line.strip().split(':')[1].split(',')
            confidences = [c for c in confstr[:-1]]
    allconfs = [int(c) for c in confidences]
    avgconf = sum(allconfs)/len(allconfs)
    a = 0.06444138 
    b = 0.39457494364395523
    return a * avgconf + b

fname = sys.argv[1]
print('JPred estimate is {:.03f}'.format(getestimate(fname)))



