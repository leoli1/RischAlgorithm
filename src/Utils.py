'''
Created on 28.09.2018

@author: Leonard
'''
from __future__ import division
import Number
import time


global log_algorithm
log_algorithm = True

numbers = [int,float,complex,Number.Rational,Number.SqrRootPolynomial]

eps = 10**(-5) # 

def Log(msg):
    if log_algorithm:
        print(msg)
def isNumber(x):
    return type(x) in numbers

def isInt(x):
    return (type(x)==int or type(x)==float or type(x)==Number.Rational) and (int(x)==x)

def numberIsZero(x):
    if type(x)==Number.SqrRootPolynomial:
        return x==0
    elif type(x)==Number.Rational:
        return x==0
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
    
def sign(x):
    return 0 if x==0 else x/abs(x)
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
    return numberIsZero(x+(-num))

def timeMethod(method):
    a = time.time()
    out = method()
    b = time.time()
    #print("{}: {}".format(method,b-a))
    return out

def addListToList(l1, l2):
    """
    l1 = [a1,a2,...,an]
    l1 = [b1,b2,...,bn]
    returns [a1+b1,a2+b2,...,an+b2]
    """
    #if len(l1)!=len(l2):
    #    raise Exception()
    
    nl = []
    for i in range(max(len(l1),len(l2))):
        nl.append(getElFromList(l1,i)+getElFromList(l2,i))
    return nl
def mulObjectToListElements(el, l):
    """
    l = [a1,a2,...,an]
    returns [a1*el,a2*el,...,an*el]
    """
    nl = []
    for e in l:
        nl.append(e*el)
    return nl
def mulListWithList(l1, l2): # cauchy product
    newDeg = len(l1)+len(l2)-2
    l = [0]*(newDeg+1)
    for i in range(newDeg+1):
        s = 0
        for j in range(i+1):
            a = getElFromList(l1, j)
            b = getElFromList(l2, i-j)
            s += a*b
        l[i] = s
    return l
def getElFromList(l,i):
    if i>=len(l):
        return 0
    return l[i]
def listStripZeros(l):
    newDeg = 0
    for i in range(len(l)-1,-1,-1):
        if not objEqualsNumber(l[i], 0):
            newDeg = i
            break
    return l[0:newDeg+1]
def listIsZero(l):
    for e in l:
        if not objEqualsNumber(e, 0):
            return False
    return True