'''
Created on 28.09.2018

@author: Leonard
'''
numbers = [int,float,complex]

eps = 10**(-7) # 

def isNumber(x):
    return type(x) in numbers

def numberIsZero(x):
    return x>-eps and x<eps