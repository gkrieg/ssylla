import sys
import os
import math
import traceback
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
#from sklearn.preprocessing import PolynomialFeatures
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
matplotlib.rcParams['axes.linewidth'] = 1.8
font = {'weight': 'semibold', 'size' : 14, 'family' : 'sans-serif'}
matplotlib.rc('font', **font)

def getconfidence(f):
    confs = []
    for line in f[2:]:
        s = line.split()
        if line[0] != '#' and len(s) == 6:
            alpha = float(s[3])
            beta = float(s[4])
            coil = float(s[5])
            l = [alpha,beta,coil]
            m = max(l)
            l.remove(m)
            nm = max(l)
            conf = 10 * (m - nm)
            
            confs.append(conf)
        elif line[0] != '#' and len(s) == 11:
            l = [float(s[a]) for a in range(3,11)]
            m = max(l)
            l.remove(m)
            nm = max(l)
            conf = 10 * (m - nm)
            confs.append(conf)
    return confs

def getq3us(test,true):
    if abs(len(test) - len(true)) == 1 and 'L' in test:
        test = test + 'L'
    elif abs(len(test) - len(true)) == 1 and 'C' in test:
        test = test + 'C'
    if len(test) == len(true):
        count = 0
        correct = 0
        lets = []
        for i,let in enumerate(test):
            count += 1
            lets.append(let)
            if let == true[i]:
                correct += 1
        return correct / count,''.join(lets)
    else:
        return -1,''

def plotconfq3(q3s,confs,weights):
    print('finallens: {} {} {}'.format(len(q3s),len(confs),len(weights)))
    halfarr = np.asarray([confs[i] for i in range(len(confs)) if q3s[i] > 0]).reshape(-1,1)
    q3arr = np.asarray([q3s[i] for i in range(len(q3s)) if q3s[i] > 0])
    weightsarr = np.asarray([weights[i] for i in range(len(q3s)) if q3s[i] > 0])
    reg = LinearRegression().fit(halfarr,q3arr,weightsarr)
    plt.scatter(confs,q3s)
    x = np.asarray([a for a in range(5,11)]).reshape(-1,1)
    y = reg.predict(x)
    plt.plot(x,y)
    plt.savefig('{}/histograms/psipred'.format(os.environ['LPGEN']))
    return reg

def geterrors(q3s,confs,reg):
    #TODO: once we get the estimator, we can actually check the error
    confarr = np.asarray(confs).reshape(-1,1)
    estimates = reg.predict(confarr)
    errors = []
    for i,est in enumerate(estimates):
        if q3s[i] > 0:
            errors.append(abs(est-q3s[i]))
    return errors

def getconfidences(prefix,r,returnarray=False,eightstate=False):
    retconfs = []
    confarrs = []
    for i in range(r):
        if r == 49:
            try:
                confs = getconfidence(open('{}{}.ss2'.format(prefix,i+1),'r').readlines())
            except:
                confs = [0]
        elif r == 292:
            try:
                confs = getconfidence(open('{}{}{}.ss2'.format(prefix,i,'6short')).readlines())
            except:
                confs = [0]
        else:
            try:
                confs = getconfidence(open('{}{}.ss2'.format(prefix,i),'r').readlines())
            except:
                confs = [0]
        if len(confs) == 0:
            print(i)
        avgconf = sum(confs)/len(confs)
        confarrs.append(confs)
        retconfs.append(avgconf)
    if returnarray == False:
        return retconfs
    else:
        return retconfs,confarrs

def getq3s(predictions,trues):
    q3s = []
    for i in range(len(predictions)):
        q3 = getq3us(predictions[i],trues[i*2+1])
        if q3 < 0:
            print(len(predictions[i]),len(trues[i*2+1]),predictions[i],trues[i*2+1])
        q3s.append(q3[0])
    return q3s

def getlens(preds):
    lens = []
    for pred in preds:
        lens.append(len(pred))
    return lens

def outputests(q3s,confs,outf,reg,lens,preds,confidencearrays):
    confarr = np.asarray(confs).reshape(-1,1)
    estimates = reg.predict(confarr)
    prettyconfs = [['{0:.2f}'.format(c) for c in conf] for conf in confidencearrays]

    for i,q in enumerate(q3s):
        if lens[i] == 1 or len(prettyconfs[i]) < 10:
            outf.write('#missing\n')
        else:
            outf.write('{} {} {} {} {}\n'.format(lens[i],q,estimates[i],preds[i].strip(),' '.join(prettyconfs[i])))


if __name__ == '__main__':
    #a test set consists of predictions,truestates, prefix of state probabilities
    testsets = [('casp10.results','{}/cleantesting/cleantestingCASP10.seqs'.format(os.environ['LPGEN']),'caspresults/casp10seq',128),('casp11.results','{}/cleantesting/cleantestingCASP11.seqs'.format(os.environ['LPGEN']),'caspresults/casp11seq',105),('casp12.results','{}/cleantesting/cleantestingCASP12.seqs'.format(os.environ['LPGEN']),'caspresults/casp12seq',55),('casp13.results','{}/cleantesting/cleantestingCASP13.seqs'.format(os.environ['LPGEN']),'caspresults/casp13/PSIPREDcasp13/CASP13',49)
            ,('6short.results','{}/cleantesting/cleantesting6short.seqs'.format(os.environ['LPGEN']),'6short/',292) 
            ,('yearly2014.results','{}/cleantesting/cleantestingyearly2014.seqs'.format(os.environ['LPGEN']),'yearly2014',100)
            ,('yearly2015.results','{}/cleantesting/cleantestingyearly2015.seqs'.format(os.environ['LPGEN']),'yearly2015',100)
            ,('yearly2016.results','{}/cleantesting/cleantestingyearly2016.seqs'.format(os.environ['LPGEN']),'yearly2016',100)
            ,('yearly2017.results','{}/cleantesting/cleantestingyearly2017.seqs'.format(os.environ['LPGEN']),'yearly2017',100)
            ,('yearly2018.results','{}/cleantesting/cleantestingyearly2018.seqs'.format(os.environ['LPGEN']),'yearly2018',100)
    #testsets = [('../caspsDeepCNF/3resultscasp10.txt','cleantestingCASP10.seqs','../caspsDeepCNF/casp10/seq',127),('resultscasp11.txt','{}/cleantesting/cleantestingCASP11.seqs'.format(os.environ['LPGEN']),'../caspsDeepCNF/casp11/seq',105),('resultscasp12.txt','{}/cleantesting/cleantestingCASP12.seqs'.format(os.environ['LPGEN']),'../caspsDeepCNF/casp12/seq',54),('resultscasp13.txt','{}/cleantesting/cleantestingCASP13.seqs'.format(os.environ['LPGEN']),'../caspsDeepCNF/casp13/seq',49)
            #,('../6MONTHS/3resultscasp10.txt','$LPGEN/cleantesting/cleantestingCASP10.seqs','casp10/casp10seq')
            ]

    q3s = []
    confidences = []
    lens = []
    preds = []
    confidencearrays = []
    weights = []
    for predictions, trues, prefix, r in testsets:
        predf = open(predictions,'r').readlines()
        preds += [p for p in predf]
        truef = open(trues,'r').readlines()
        lens += getlens(predf)
        print(len(predf),len(truef))
        q3s += getq3s(predf,truef)
        confs,confarr = getconfidences(prefix,r,returnarray=True,eightstate=True)
        confidences += confs
        confidencearrays += confarr
        if 'CASP' in trues:
            weight = .5/4*2/len(predf)
        else:
            weight = .5/6/len(predf)*2
        #confidences += getconfidences(prefix,r)
        weights += [weight for w in predf]
        print('q3s: {} confidences {}'.format(len(q3s),len(confidences)))
        print(sum(q3s)/len(q3s))

    reg = plotconfq3(q3s,confidences,weights)
    errs = geterrors(q3s,confidences,reg)
    print('error is {}'.format(sum(errs)/len(errs)))
    outputests(q3s,confidences,open('psipredestimates.txt','w'),reg,lens,preds,confidencearrays)
        

