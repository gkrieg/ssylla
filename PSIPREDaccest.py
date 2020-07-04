import sys

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


if __name__ == '__main__':
    confs = getconfidence(open(sys.argv[1],'r').readlines())
    avgconf = sum(confs)/len(confs)
    a = 0.06684236
    b = 0.32518820361974776
    estimate = a * avgconf + b
    print('PSIPRED estimated accuracy is {:.03f}'.format(estimate))

        

