# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 14:27:55 2019

@author: Dan

region.py scores strands based on their region lengths
"""

# shortest region on strand
# returns integer count
def calc_short_region(strand):
    reglen = 0 # region length
    minreglen = len(strand) # tracks smallest seen region length
    curr = None # current ring
    for n in strand:
        this = n.up().up() # ring of nucleotide
        if curr == None: # if not set, for 1st base in strand
            curr = this
            reglen+=1
        elif curr == this: # if still on the same ring
            reglen+=1
        elif curr != this: # if the nucl moved to a new ring
            minreglen = min(minreglen,reglen) # terminate the region
            reglen = 1 # restart the counter
            curr = this # set the new current ring
    # need one more to catch when strand terminates
    minreglen = min(minreglen,reglen) # terminate the region
    
    return minreglen

# strand has region less than length
# returns binary
def has_short_region(strand,length):
    if calc_short_region(strand)<length:
        return True
    else:
        return False