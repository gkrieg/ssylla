import sys

def get3estimate(fname):
    inf = open(fname).readlines()
    confidences = []
    for line in inf:
        if '#' not in line:
            probstr = line.strip().split(',')[:-1]
            stateprobabilities = [float(c) for c in probstr]
            maxprob = max(stateprobabilities)
            stateprobabilities.remove(maxprob)
            secondmaxprob = max(stateprobabilities)
            confidence = maxprob - secondmaxprob
            confidences.append(confidence)
    avgconf = sum(confidences)/len(confidences)
    a = 0.06354173
    b = 0.36536894983769763
    return a * avgconf + b
  
def get8estimate(fname):
    inf = open(fname).readlines()
    confidences = []
    for line in inf:
        if '#' not in line:
            probstr = line.strip().split(',')[:-1]
            stateprobabilities = [float(c) for c in probstr]
            maxprob = max(stateprobabilities)
            stateprobabilities.remove(maxprob)
            secondmaxprob = max(stateprobabilities)
            confidence = maxprob - secondmaxprob
            confidences.append(confidence)
    avgconf = sum(confidences)/len(confidences)
    a = 0.07204353
    b = 0.2806027465372139
    return a * avgconf + b

fname = sys.argv[1]
numstructurestates = int(sys.argv[2])
if numstructurestates == 3:
    print('MUFOLD estimate is {:.03f}'.format(get3estimate(fname)))
elif numstructurestates == 8:
    print('MUFOLD estimate is {:.03f}'.format(get8estimate(fname)))


