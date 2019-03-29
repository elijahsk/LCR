#!/usr/bin/python2.7

import sys, random, math;
from operator import itemgetter

### Generates a random directed graph with edge labels using preferential attachment
### The labels can follow either a uniform or normal distribution

def exitMain():
    sys.exit(1);

if(len(sys.argv) != 4):
    print("usage python2.7 " + sys.argv[0] + " unlabeled_file_path L {uni|norm|exp}");
    exitMain();

try:
    filePath = sys.argv[1];
    L = int(sys.argv[2]);
    dist = sys.argv[3];
    name = filePath + "_labeled";
except:
    print("second param need to be integers");
    exitMain();


infile = open(filePath, 'r');
f = open(name, 'w');


random.seed();
mean = int(math.floor(L/2));
sd = int( max(1, math.floor(L/4)) );

dgraph = [];
line = infile.readline();
while line:
    nodes = line.split(' ');
    node1 = nodes[0];
    node2 = nodes[1]; 

    if dist == "norm":
        label = random.normalvariate(mean, sd);
    if dist == "uni":
        label = random.uniform(0,L);
    if dist == "exp":
        label = random.expovariate( 1.0 / L/1.7 )
        

    label = max(0,label);
    label = min(L-1,label);
    label = int(math.floor(label));

    #label = 2**label;
    #print(label);
       
    dgraph.append((node1, node2, label));

    line = infile.readline();
infile.close();

dgraph = sorted(dgraph,key=lambda x: x[0], reverse=False);

for triple in dgraph:
    (s,t,l) = triple;
    base = 0;
    ss = str(s) + " " + str(t) + " " + str(l);
    f.write(ss + "\n");

f.close();
