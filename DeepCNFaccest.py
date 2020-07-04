import sys

def getconfidence(f):
    confs = []
    for line in f:
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

if __name__ == '__main__':
    f = open(sys.argv[1],'r').readlines()
    confs = getconfidence(f)
    avgconf = sum(confs)/len(confs)
    if sys.argv[2] == '3state':
        a = 0.04219216 
        b = 0.557812480435692
    elif sys.argv[2] == '8state':
        a = 0.07449779
        b = 0.2470244144005686
    estimate = a * avgconf + b
    print('DeepCNF estimate is {:.03f}'.format(estimate))
        

