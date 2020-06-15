import sys

def getestimatedist(f,prefix,returnarray=False,confidences=None):
    #Fit parameters for the piecewise-linear fit that make up the percent identity estimator
    a = 2.5 
    b = -15 
    c = 0.2288796586750804 
    d = 71.80858692549229
    e = 39

    #getting the per-residue estimates from the percent identity file
    accs = []
    for i,val in enumerate(f):
        val = float(val.strip())
        boxacc = getdistacc(a,b,c,d,val,e)
        accs.append(boxacc)

    #getting the piecewise-linear fit of the confidences
    conf = sum(confidences)/len(confidences)
    conf = processconfidences([conf],prefix=prefix)[0]
    confs = processconfidences(confidences,prefix=prefix)
    boxacc = sum(accs)/len(accs)

    #combining the two into a final estimate and a per-residue estimate
    if pi == True:
        if prefix == '8state':
            x = 0.70530552 
            y = 0.5342837
            z = -0.20529427223651175
        else:
            x = 0.67104363 
            y = 0.39915731
            z = -0.05542468172164372
    acc = conf * x + boxacc * y + z
    accs = [confs[i] * x + accs[i] * y + z for i in range(len(confidences))]
    return acc,accs

def func(x,a,b):
    return a * x + b

def getdistacc(a,b,c,d,val,e):
    return func(val*100,a,b)/100 if val * 100 < e else func(val*100,c,d)/100

def getnnessyconfidences(f):
    confs = []
    for line in f:
        probs = [float(s) for s in line.split()]
        mprob = max(probs)
        probs.remove(mprob)
        m2prob = max(probs)
        conf = (mprob - m2prob) * 10
        confs.append(conf)
    return confs

def processconfidences(confidences,prefix='3state'):
    retconfs = []
    for conf in confidences:
        if conf < 4:
            if prefix == '8state':
                retconf = conf * 0.20401978 +  0.07444265879675022
            else:
                retconf = conf * 0.1972279 + 0.07767501650766317
            retconfs.append(retconf)
        else:
            if prefix == '8state':
                retconf = conf * 0.02307881 + 0.7433726706040007
            else:
                retconf = conf * 0.01845465 +  0.7834331948971939
            retconfs.append(retconf)
    return retconfs
            
def getnnessyestimate():
    prefix = sys.argv[3]
    confs = getnnessyconfidences(open(sys.argv[1],'r').readlines())
    perif = open(sys.argv[2],'r').readlines()
    estimate,estimatearr = getestimatedist(perif,prefix,returnarray=True,confidences=confs)
    return estimate

if __name__ == '__main__':
    print(getnnessyestimate())

