'''
Created on 28.09.2018

@author: Leonard
'''
import Number

numbers = [int,float,complex,Number.Rational,Number.SqrRootPolynomial]

eps = 10**(-5) # 

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
