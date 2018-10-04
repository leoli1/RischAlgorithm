'''
Created on 28.09.2018

@author: Leonard
'''
numbers = [int,float,complex]

eps = 10**(-5) # 

def isNumber(x):
    return type(x) in numbers
def isInt(x):
    return (type(x)==int or type(x)==float) and (int(x)==x)

def numberIsZero(x):
    return abs(x)<eps

def argmax(l):
    r = range(len(l))
    return max(r, key=lambda i: l[i])
def isPoly(x):
    try:
        x.degree
        return True
    except AttributeError:
        return False
    
def objEqualsNumber(obj, num):
    if obj==None:
        return False
    x = None
    if isNumber(obj):
        x = obj
    else:
        y = obj.getConstant()
        if y==None:
            return False
        x = y
    return numberIsZero(x-num)
