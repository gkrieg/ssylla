import sys

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

def getpercentidentitymatch(indf):
    currprot = ''
    currlist = []
    for line in indf:
        if '>' in line:
            currprot = line[1:].strip()
        elif len(line.split()) == 11:
            currlist.append((float(line.split()[3]),float(line.split()[0]),float(line.split()[4])))
    return [sum([x[0] for x in currlist])/len(currlist)]

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
     
    newcol,newcolarr,maxprob = getprob(open(sys.argv[1],'r').readlines(),returnarray=True,maxprob=True)
    sub = transformsubs(getpercentidentitymatch(open(sys.argv[2],'r').readlines()))[0]
    if sys.argv[3] == '8state':
        a = 0.97222608 
        b = -0.03719895  
        c = 1.01279835 
        d = -0.41289329398929386
    elif sys.argv[3] == '3state':
        a = 0.51497103 
        b = -0.14905368  
        c = 3.01569142 
        d = -0.9629875662193697
    estimate = a * sub + b * newcol[0] + c * maxprob[0] + d
    print('SSpro estimate is {:.03f}'.format(estimate))
