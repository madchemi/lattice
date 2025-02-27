import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import math 
from math import *
import random


def read_init(filename):
    """
    init.inp을 읽고 해당하는 값들을 반환 
    
    init.inp structure
     ----  Degree of Polymerization
     ----  time 
     ----  number of bins 
     ----  maximum distance of x 
    
    """
    f = open('init.inp', mode = 'r')
    lines = f.readlines()
    text = []
    print(lines)
    for line in lines:
        line = line.strip()
        text.append(line)
    f.close()
    ndeg, ntime, nbin, amax = text[0], text[1], text[2],text[3]

    return ndeg, ntime, nbin, amax

def make_data(data):
    """
    read_init -> tuple 형태의 데이터를 반환 
    따라서, 데이터를 받은 다음에 tuple에서 해당하는 변수를 가져와야함. 
        init.inp structure
     ----  Degree of Polymerization
     ----  time 
     ----  number of bins 
     ----  maximum distance of x 
    변수명은 위와 동일하게 

    """
    ndeg, ntime, nbin, amax = int(data[0]), int(data[1]),int(data[2]),int(data[3])
    dist = [0]*nbin 
    bin = float(amax/nbin)
    total_ete  = 0.0
    
    for i in range(ntime):
        x0,y0,z0 = 0.0,0.0,0.0
        a_ete = 0.0
        for j in range(ndeg):
            a = random.randint(1,6)
            if a == 1: 
                x0  = x0 + 1.0
            elif a == 2:
                x0 = x0 -1.0
            elif a == 3:
                y0 = y0 + 1.0
            elif a == 4:
                y0 = y0 -1.0
            elif a ==5:
                z0 = z0 + 1.0
            elif a == 6:
                z0 = z0 - 1.0
        a_ete = sqrt(x0*x0 + y0*y0 + z0*z0)
        total_ete  = total_ete + a_ete**2.0
        inn = int(a_ete/bin) + 1
        if (inn<nbin):
            dist[inn] = dist[inn] + 1
    aver_ete = total_ete/float(ntime)
    print(f'DOF: {ndeg}, ete: {aver_ete}') 


    f = open(f'./data/dist.out.{ndeg}', mode = 'w')
    for i in range(nbin):
        x = bin/2.0 + float((i-1))*bin 
        x1 = x - bin/2.0
        x2 = x + bin/2.0 
        dnorm = 1.0 #4.0*pi(x2**3.0 - x1**3.0)/3.0 
        y = dist[i]/dnorm/float(ntime)/bin

        f.write(f'{x} {y} \n')
    f.close()

    f2 = open(f'./data/ref.out.{ndeg}', 'w')
    npoint = 500
    bin_r = float(amax)/float(npoint)
    for i in range(npoint):
        x = float(i-1)*bin_r
        y = 4.0*pi*(3.0/(2.0*pi*float(ndeg)))**(3.0/2.0)*exp(-3.0*x*x/2.0/float(ndeg))*x*x
        f2.write(f'{x} {y} \n')
    f2.close()
asd = read_init('init.inp')    
make_data(asd)
             

