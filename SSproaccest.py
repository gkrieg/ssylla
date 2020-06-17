import sys
import os

import math
import traceback
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
#from sklearn.preprocessing import PolynomialFeatures
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
matplotlib.rcParams['axes.linewidth'] = 1.8
font = {'weight': 'semibold', 'size' : 14, 'family' : 'sans-serif'}
matplotlib.rc('font', **font)
import numpy as np
from scipy import stats


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

def clean(a):
    arr = []
    for l in a:
        if l in 'GHIBELTS':
            arr.append(l)
        else:
            arr.append('L')
    return ''.join(arr)
def clean3(a):
    arr = []
    for l in a:
        if l == 'H':
            arr.append('A')
        elif l == 'E':
            arr.append('B')
        else:
            arr.append('C')
    return ''.join(arr)

def gethijackedscore(hf):
    total = int(hf[0][1:])
    subbed = int(hf[2])
    return subbed / total

def getentropy(entf):

    ents = []
    ent = []
    for line in entf[1:]:
        l = line.split()
        if len(l) == 1:
            ents.append(sum(ent)/len(ent))
            ent = []
        else:
            #they are 5,6,7
            #process end of protein
            s = l[5:]
            entropy = -sum([float(t) * math.log(float(t)) for t in s if float(t) != 0])
            ent.append(entropy)
    ents.append(sum(ent)/len(ent))
    return ents

def getprob(entf,returnarray=False,maxprob=False):

    ents = []
    entarray = []
    ent = []
    topss = []
    tops = []
    for line in entf[1:]:
        l = line.split()
        if len(l) == 1:
            ents.append(sum(ent)/len(ent))
            entarray.append(ent)
            ent = []
            topss.append(sum(tops)/len(tops))
            tops = []
        else:
            #they are 5,6,7
            #process end of protein
            s = [float(a) for a in l[5:]]
            top = max(s)
            tops.append(top)
            s.remove(top)
            second = max(s)
            ent.append((top - second) * 10)
    ents.append(sum(ent)/len(ent))
    entarray.append(ent)
    topss.append(sum(tops)/len(tops))
    if returnarray == False:
        return ents
    else:
        if maxprob == False:
            return ents,entarray
        else:
            return ents,entarray,topss

def getpercentidentitymatch(dirr):
    #given the directory, we need to go through all the index files
    proteins = {}
    for f in os.listdir(dirr):
        if 'index' in f:
            indf = open('{}/{}'.format(dirr,f),'r').readlines()
            #now we just process the files
            currprot = ''
            for line in indf:
                if '>' in line:
                    currprot = line[1:].strip()
                    proteins[currprot] = []
                elif len(line.split()) == 11:
                    proteins[currprot].append((float(line.split()[3]),float(line.split()[0]),float(line.split()[4])))
    #TODO:Need to now take the proteins dictionary and really split it up so that each protein is in order and we can return any of them so that we can actually plot them finally.
    plist = []
    for i in range(len(proteins.keys())):
        key = 'p{}'.format(i)
        avglist = []
        for prot in proteins[key]:
            avglist.append(prot[0])
            #avglist.append(prot[1])

        if len(avglist) > 0:
            plist.append(sum(avglist)/len(avglist))
        else:
            plist.append(0)

    return plist

def transformsubs(subs):

    retsubs = []
    a = 0.00814905 
    b = 0.39277664467985623

    c = 0.00266278 
    d = 0.6507551242336596
    for sub in subs:
        if sub > 50:
            retsub = c * sub + d
        elif sub < 5:
            retsub = .6
        else:
            retsub = a * sub + b
        retsubs.append(retsub)
    return retsubs




if __name__ == "__main__":
    import os
    testsets = [('ocasp102.ss','{}/cleantesting/cleantestingCASP10.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200330-162500-300857245041/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200330-162404-817303602682/sspro.out')
            ,('casp112.ss','{}/cleantesting/cleantestingCASP11.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200406-160706-670936510996/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200406-160628-201220748113/sspro.out')
            ,('casp122.ss','{}/cleantesting/cleantestingCASP12.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200406-165233-285414849506/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200406-165211-089901953102/sspro.out')
            ,('casp132.ss','{}/cleantesting/cleantestingCASP13.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200406-173758-858640086935/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200406-173718-583940793858/sspro.out')
            ,('6short.ss','{}/cleantesting/cleantesting6short.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200508-161713-855151926084/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200508-161433-662844601004/sspro.out')
            ,('oyearly2014.ss','{}/cleantesting/cleantestingyearly2014.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-112042-429503635601/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200520-111947-405358173983/sspro.out')
            ,('oyearly2015.ss','{}/cleantesting/cleantestingyearly2015.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-113249-263805504055/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200520-113148-481019451333/sspro.out')
            ,('oyearly2016.ss','{}/cleantesting/cleantestingyearly2016.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-132518-240452903421/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200520-132413-580096768132/sspro.out')
            ,('oyearly2017.ss','{}/cleantesting/cleantestingyearly2017.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-133427-080019286932/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200520-133322-822651519538/sspro.out')
            ,('oyearly2018.ss','{}/cleantesting/cleantestingyearly2018.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-145241-098144359053/'.format(os.environ['SSPRO']),'../pkg/SSpro_5.2/tmp/20200520-145152-864507009373/sspro.out')]

    testsets = [('ocasp102.ss8','{}/cleantesting/cleantesting8stateCASP10.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200330-162500-300857245041/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200504-172838-471937589264/sspro8.out'),('casp112.ss8','{}/cleantesting/cleantesting8stateCASP11.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200406-160706-670936510996/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200504-172339-900717801513/sspro8.out'),('casp122.ss8','{}/cleantesting/cleantesting8stateCASP12.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200406-165233-285414849506/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200504-182457-602502753125/sspro8.out'),('casp132.ss8','{}/cleantesting/cleantesting8stateCASP13.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200406-173758-858640086935/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200504-184300-745095509203/sspro8.out'),('6short.ss8','{}/cleantesting/cleantesting8state6short.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200508-161713-855151926084/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200508-161433-281607839377/sspro8.out')
            ,('oyearly2014.ss8','{}/cleantesting/cleantesting8stateyearly2014.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-112042-429503635601/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200520-111947-333854095177/sspro8.out')
            ,('oyearly2015.ss8','{}/cleantesting/cleantesting8stateyearly2015.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-113249-263805504055/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200520-113148-552042418258/sspro8.out')
            ,('oyearly2016.ss8','{}/cleantesting/cleantesting8stateyearly2016.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-132518-240452903421/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200520-132413-713057660596/sspro8.out')
            ,('oyearly2017.ss8','{}/cleantesting/cleantesting8stateyearly2017.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-133427-080019286932/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200520-133322-158988551335/sspro8.out')
            ,('oyearly2018.ss8','{}/cleantesting/cleantesting8stateyearly2018.seqs'.format(os.environ['LPGEN']),'{}/pkg/HOMOLpro_1.1/tmp/20200520-145241-098144359053/'.format(os.environ['SSPRO']),'../pkg/SSpro8_5.2/tmp/20200520-145152-326374908076/sspro8.out')]
      
     
    col = []
    colinput = []
    colinputarray = []
    q3sa = []
    subs= []
    lens = []
    predseqs = []
    confs = []
    maxprobs = []
    weights = []
    for inp,trues,subin,colins in testsets:
        inpreds = open(inp,'r').readlines()
        trueseqs = open(trues,'r').readlines()
        inpredslen = len(inpreds)//2
        q3s = []
        if 'CASP' in trues:
            weight = .5 / 4 / len(trueseqs) * 2
        else:
            weight = .5 / 6 / len(trueseqs)*2
        for i in range(inpredslen):
            inpred = clean(inpreds[i*2+1])
            #inpred = clean3(inpreds[i*2+1])
            predseqs.append(inpred)
            lens.append(len(inpred))
            q3 = getq3us(inpred,trueseqs[i*2+1])
            q3s.append(q3)
            weights.append(weight)

            #get the hijacked numbers of how many residues were scored by homology
            #hijackedf = open('p{}.hijacked'.format(i),'r').readlines()
            #hinum = gethijackedscore(hijackedf)
            #subs.append(hinum)

        #get entropy
        #subs = getentropy(open(sys.argv[3],'r').readlines())
        cm = plt.get_cmap("RdYlGn")
        newcol,newcolarr,maxprob = getprob(open(colins,'r').readlines(),returnarray=True,maxprob=True)
        maxprobs = maxprobs + maxprob
        colinput = colinput + newcol
        colinputarray = colinputarray + newcolarr
        #colinput = colinput + getprob(open(colins,'r').readlines()) 
        col = col + [cm(colin / 10) for colin in colinput]
        subs = subs + getpercentidentitymatch(subin)
        print(len(subs),len(q3s),len(maxprobs))
        print(subs[0],colinput[0])

        q3sa = q3sa + [q3[0] for q3 in q3s]
    numplots = 5
    #for i in range(numplots + 1):
    i = 2
    plt.figure(i)
    subs = transformsubs(subs)
    features = np.asarray([[subs[j],colinput[j],maxprobs[j]] for j in range(len(subs))])
    #features = np.asarray([[subs[j],colinput[j]] for j in range(len(subs))])
    #halfnhalf = [.01/numplots * i * subs[j] + colinput[j] * .1/numplots * (numplots - i) for j in range(len(subs))]
    #halfnhalf = np.asarray([subs[j] for j in range(len(subs)) if subs[j] > 20 and subs[j] < 50]).reshape(-1,1)
    halfnhalf = np.asarray([colinput[j] for j in range(len(colinput))]).reshape(-1,1)
    #halfnhalf = np.asarray([subs[j] for j in range(len(subs))]).reshape(-1,1)
    #q3sb = np.asarray([q3sa[j] for j in range(len(subs)) if subs[j] > 20 and subs[j] < 50])
    #halfnhalf2 = np.asarray([subs[j] for j in range(len(subs)) if subs[j] > 50]).reshape(-1,1)
    #q3sc = np.asarray([q3sa[j] for j in range(len(subs)) if subs[j] > 50])
    #hhmin = min(halfnhalf)
    #hhmax = max(halfnhalf) - hhmin
    #confs = [(10/hhmax) * (h - hhmin) for h in halfnhalf]
    #plt.scatter(halfnhalf,q3sa,c=col)
    plt.scatter(halfnhalf,q3sa)
    #plt.scatter(halfnhalf2,q3sa)

    reg1 = LinearRegression().fit(halfnhalf,q3sa,weights)
    print(reg1.coef_,reg1.intercept_)
    reg1.predict(halfnhalf)
    x = np.asarray([a for a in range(10)]).reshape(-1,1)
    y = reg1.predict(x)
    plt.plot(x,y)
    #reg2 = LinearRegression().fit(halfnhalf2,q3sc)
    #print(reg2.coef_,reg2.intercept_)
    #reg2.predict(halfnhalf2)
    #x = np.asarray([a for a in range(50,100)]).reshape(-1,1)
    #y = reg2.predict(x)
    #plt.plot(x,y)

    halfarr = np.asarray(halfnhalf)
    halfarr = halfarr.reshape(-1,1)
    q3arr = np.asarray(q3sa)
    #reg = LinearRegression().fit(halfarr,q3arr)
    reg = LinearRegression().fit(features,q3arr)
    outf = open('ssproestimates.txt','w')
    #preds = reg.predict(halfarr)
    preds = reg.predict(features)
    errors=[]
    for i in range(len(q3arr)):
        #TODO: Need to get the and confidences and estimates for 8 state
        q3 = q3arr[i]
        pred = preds[i]
        l = lens[i]
        predseq = predseqs[i]
        conf = ['{0:.2f}'.format(c) for c in colinputarray[i]]
        outf.write('{} {} {} {} {}\n'.format(l,q3,pred,predseq,' '.join(conf)))
        error = abs(q3-pred)
        errors.append(error)


    #plt.plot(x,y) 
    plt.savefig('{}/histograms/ssprohybrid'.format(os.environ['LPGEN']))
    print('Average error: {}'.format(sum(errors)/len(errors)))
