import sys
import os
from sklearn.linear_model import LinearRegression
import numpy as np
import traceback
outf = open('jpredestimates.txt','w')

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



def convert(a):
    r = []
    for l in a:
        if l == 'H':
            r.append('A')
        elif l == 'E':
            r.append('B')
        elif l != '\n':
            r.append('C')
    return ''.join(r)

def getestimates(allq3s,confidences,weights):
    allconfstr = [[int(c) for c in con] for con in confidences]
    confstr = [[int(c) for c in con] for con in confidences if len(con) > 5]
    avgconfs = np.asarray([sum(conf)/len(conf) for conf in confstr]).reshape(-1,1)
    allavgconfs = np.asarray([sum(conf)/len(conf) for conf in allconfstr]).reshape(-1,1)
    q3s = np.asarray([allq3s[i] for i in range(len(allq3s)) if len(confidences[i]) > 5])
    weightsarr = np.asarray([weights[i] for i in range(len(allq3s)) if len(confidences[i]) > 5])
    reg = LinearRegression().fit(avgconfs,q3s,weightsarr)
    pred = reg.predict(allavgconfs)
    return pred




testsets = [('casp10',128),('casp11',105),('casp12',55),('casp13',49),('6MONTHSnew',309)
        #,('yearly2014',100)
        ,('yearly2015',100),('yearly2016',100),('yearly2017',100),('yearly2018',100)]
predictions = []
q3s = []
confidences = []
estimates = []
weights = []
for testset,r in testsets:
    if testset == '6MONTHSnew' or 'yearly' in testset:
        seqlines = open('{}/cleantesting/cleantesting{}.seqs'.format(os.environ['LPGEN'],testset),'r').readlines()
    else:
        seqlines = open('{}/cleantesting/cleantesting{}.seqs'.format(os.environ['LPGEN'],testset.upper()),'r').readlines()
    offset = 0
    for i in range(r):
        try:
            w = .5/len(seqlines)/2/4*2 if 'casp' in testset else .5/len(seqlines)/2/5*2
            weights.append(w)
            if testset == '6MONTHSnew':
                inf = open('{}_output/seq{}.jnet'.format('6MONTHS',i),'r').readlines()
            else:
                inf = open('{}_output/seq{}.jnet'.format(testset,i),'r').readlines()
            skipconfidence = False
            for line in inf:
                if 'jnetpred' in line:
                    if testset == '6MONTHSnew':
                        seqline = seqlines[(i-offset)*2+1]
                    else:
                        seqline = seqlines[i*2+1]
                    l = line.split(':')[1].split(',')
                    r = convert(l)
                    if len(seqline) != len(r) + 1:
                        offset += 1
                        skipconfidence = True
                    else:
                        predictions.append(r)
                        q3,s = getq3us(r,seqline)
                        q3s.append(q3)
                elif 'JNETCONF' in line and skipconfidence == False:
                    confstr = line.strip().split(':')[1].split(',')
                    confs = [c for c in confstr[:-1]]
                    confidences.append(confs)
        except:
            print('seq{} failed'.format(i))
            q3s.append(0)
            confidences.append(['0'])
            predictions.append('')
            traceback.print_exc()
    print(len(q3s),len(confidences),len(weights))
estimates = getestimates(q3s,confidences,weights)
print(len(q3s),len(confidences),len(predictions),len(estimates))
for i in range(len(q3s)):
    if q3s[i] == 0:
        outf.write('missing\n')
    else:
        outf.write('{} {} {} {} {}\n'.format(len(predictions[i]),q3s[i],estimates[i],predictions[i],' '.join(confidences[i])))
