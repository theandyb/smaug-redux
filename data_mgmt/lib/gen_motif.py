from itertools import product

alpha = ['A', 'T', 'C', 'G']
motifs = [a+b+c+d+e+f+g for a,b,c,d,e,f,g in product(alpha, repeat=7)]
with open('motifs7.txt', 'w') as f:
    for item in motifs:
        f.write("%s\n" % item)
