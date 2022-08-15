# -*- coding: utf-8 -*-
"""
Created on Thu May 9 2019

@author: Dan

seeding.py

part of debug module

analysis of seeding regions
"""

# calculate the average of a 1-D list
def mean(li):
    return sum(li)/len(li)

class Mixin:
    # returns a single score
    def seedscore(self,seedlen):
        o = minseedlen(self,seedlen)
        return round(mean(o),4)
    
#    # returns a score based on how many seed sites each staple has
#    # this is probably crap
#    # don't use it
#    def seedperscore(self,seedlen):
#        o = numseedlen(self,seedlen)
#        return round(mean(o),4)
    
    # returns a single score for consolidated section
    def ringseedscore(self,seedlen,selection):
        o = ringminseedlen(self,seedlen,selection)
        return mean(o)

# checks how many staple strands have a minimum seeding region 'seedlen'
# returns a list of binary values
def minseedlen(origami,seedlen):
    all_staples = origami.get_all_staples()
    seed_list = []
    for s in all_staples:
        seed_list.append(has_seed(s,seedlen))
    return seed_list
    
# returns a score for a consolidated section of staples
def ringminseedlen(origami,seedlen,selection):
    selected_staples = origami.get_ring_staples(selection)
    seed_list = []
    for s in selected_staples:
        seed_list.append(has_seed(s,seedlen))
    return seed_list

## returns a score based on how many seed sites each staple has
## this is probably crap
## don't use it
#def numseedlen(origami,seedlen):
#    all_staples = origami.get_all_staples()
#    seed_list = []
#    for s in all_staples:
#        seed_list.append(num_seed(s,seedlen))
#    return seed_list
#
## checks the regions on a staple strand and counts how many nucleation seeds there are
## this is probably crap
## don't use it
#def num_seed(strand,seedlen):
#    reglen = 0 # region length
#    seednum = 0  # number of seeds
#    curr = None # current ring
#    for n in strand:
#        this = n.up().up() # ring of nucleotide
#        if curr == None: # if not set, for 1st base in strand
#            curr = this
#            reglen+=1
#        elif curr == this: # if still on the same ring
#            reglen+=1
#        elif curr != this: # if the nucl moved to a new ring
#            if reglen>=seedlen: # check if region is a seed
#                seednum+=1
#            reglen = 1 # restart the counter
#            curr = this # set the new current ring
#    if reglen>=seedlen: # check the last region
#        seednum+=1
#    return seednum

# checks a single staple for the presence of a nucleation region
def has_seed(strand,seedlen):
    reglen = 0 # region length
    curr = None # current ring
    for n in strand:
        this = n.up().up() # ring of nucleotide
        if curr == None: # if not set, for 1st base in strand
            curr = this
            reglen+=1
        elif curr == this: # if still on the same ring
            reglen+=1
        elif curr != this: # if the nucl moved to a new ring
            if reglen>=seedlen: # check if region is a seed
                return True
            reglen = 1 # restart the counter
            curr = this # set the new current ring
    if reglen>=seedlen: # check the last region
        return True
    return False