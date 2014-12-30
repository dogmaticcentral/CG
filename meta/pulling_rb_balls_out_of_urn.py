#!/usr/bin/python -u


import sys, os, math, random
from scipy import special

def beta_binomial (n1, N1, n2, N2):
    retval = 0
    retval += special.gammaln(N1+1) - special.gammaln(n1+1)
    retval += special.gammaln(N2+2) - special.gammaln(n2+1)
    retval += special.gammaln(n1+n2+1) - special.gammaln(N1+N2+2)
    retval +=  special.gammaln( N1 + N2 - n1 - n2 + 1) - special.gammaln(N1-n1+1) - special.gammaln(N2-n2+1) 
    retval = math.exp(retval)
    return retval

#########################################
def main():

    # number of balls in a bag
    N    = 450
    # Blue of which are of the color blue, the rest are red
    Blue = 70
    
    # we are pulling  balls out at random
    # the probability  of succesfully pulling the ball out is mu
    mu = 0.01

    # if we repeat the drawing experiment very many times, 
    # what is the expected number of blue balls we'll manage to pull out?
    expected = 0
    for blue  in range (Blue+1):
        for n  in range (blue, N+1): # n is the number of balls drawn
            expected += special.binom(Blue, blue)*special.binom(N-Blue, n-blue)*math.pow(mu,n)*math.pow(1-mu, N-n)*blue
    print "expected, calculated: ", expected

    expectedMIT = 0
    Red = N-Blue
    
    for  blue  in range (Blue+1):
        expectedMIT += blue*beta_binomial (blue, Blue, mu*Red, Red)
    print "expectedMIT, calculated: ", expectedMIT
   
    # now try the same thing through a simulation
    bag = []
    for i in range(N): 
        bag.append("r")

    for i in range(Blue):
        done = False
        while not done:
            rand_pos = random.randint(0,N-1)
            if bag[rand_pos] == "r":
                bag[rand_pos] = "b"
                done = True

    avg = 0.0
    # say we do 1000 experiments
    number_of_exps = 1000
    for exp in range(number_of_exps):
        count = 0
        for i  in range(N):
            if random.random() < mu:
                # we have chosen this position ("pulled out this ball")
                if bag[i] == "b": count +=1
        avg += count
    avg /= number_of_exps
    print "avg in the experiment: ", avg
   



#########################################
if __name__ == '__main__':
    main()
