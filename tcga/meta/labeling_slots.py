#!/usr/bin/python -u

#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#


import sys, os, math, random
from scipy import special
from random import randrange

#########################################
def simulation (M, Nr, Nb, l, number_of_iterations):
    
    avg_number_of_double_labeled = 0
    pval = 0.0

    if not number_of_iterations > 0:
        return  [avg_number_of_double_labeled, pval]


    for i in range(number_of_iterations):
        #####
        slots = []
        for s in range(M):
            slots.append({"r":0, "b":0})
        number_of_double_labeled = 0
        for j in range(Nr):
            random_slot = randrange(M)
            slots[random_slot]["r"] += 1
        for j in range(Nb):
            random_slot = randrange(M)
            slots[random_slot]["b"] += 1

        for s in range(M):
            if slots[s]["r"]>0  and  slots[s]["b"]>0:
                #print " %3d   %2d  %2d " %  (s, slots[s]["r"] ,  slots[s]["b"])
                number_of_double_labeled += 1

        #####
        avg_number_of_double_labeled += number_of_double_labeled
        if ( number_of_double_labeled <= l ): pval += 1.0

    ##################################
    avg_number_of_double_labeled /= float(number_of_iterations)
    pval /= float(number_of_iterations)

    return [avg_number_of_double_labeled, pval]
        
#########################################
def main():

    # we have M slots
    M = 6996
    # we are throwing in Na red balls (labels)
    # balls perhaps are better analogy because  we can have two balls/labels in the same slot
    Nr = 2907
    # for black label
    Nb = 62

    # what is the expected number of slots, <m>, that have both red and black label?
    # what is the probability of having l or less slots with both label?
    l = 13
    number_of_iterations = 10000
    [avg, prob]  = simulation (M, Nr, Nb, l, number_of_iterations)

    print " avg  %6.2f   pval:%6.3f " % (avg, prob)
    
    
   



#########################################
if __name__ == '__main__':
    main()
