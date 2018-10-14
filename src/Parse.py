'''
Created on 28.09.2018

@author: Leonard
'''
from __future__ import division

import Polynomial as Pol
import RationalFunction as Rat
import FieldExtension as FE
from Utils import *
import Number
import time


def parseField0PolyFromStr(poly_raw,var="x"):
    monomials = poly_raw.split("+")
    poly = Pol.Polynomial()
    for mon in monomials:
        if len(mon)==0:
            continue
        l = mon.split("**")
        i = complex(0,1)
        power = 0
        if (len(l)>1):
            power = eval(l[1])
        elif l[0].endswith(var):
            power = 1
        coeff_raw = l[0].strip(var).strip("*")
        coeff = 1
        if len(coeff_raw)!=0:
            coeff = eval(coeff_raw)
        poly += Pol.Monomial(power,coeff)
    return poly
def parseField0RatFromStr(rat_raw,var="x"):
    parts = rat_raw.split("/")
    if len(parts)>2:
        raise Exception()
    if len(parts)==1:
        return 
    return Rat.RationalFunction(parseField0PolyFromStr(parts[0], var), parseField0PolyFromStr(parts[1], var))

def parseExpressionFromStr(expr_str,fieldTower=None):
    a = time.time()
    if fieldTower==None:
        fieldTower = FE.BASEFIELD
    variables = [fieldTower.getFieldExtension(i).variable for i in range(fieldTower.towerHeight)]#[(fieldTower.getFieldExtension(i).variable, fieldTower.getStrippedTower(i+1)) for i in range(fieldTower.towerHeight)]

    for variable in variables:
        poly = Pol.Polynomial([Number.ZERO,Number.ONE],variable=variable)#fieldTower=tower)
        exec (variable.stringRepr + "=poly")
    i = complex(0,1)
    try:
        expr = eval(expr_str)
        if isNumber(expr):
            expr = Pol.Polynomial([expr],variable=FE.BASEVARIABLE)
       # if isNumber(expr) or expr.fieldTower.towerHeight<fieldTower.towerHeight:
       #     expr = Pol.Polynomial([expr],fieldTower=fieldTower)
            
        expr.replaceNumbersWithRationals()
        b = time.time()
        print("parsing expression took: {}s".format(b-a))
        return expr
    except (SyntaxError,NameError,ZeroDivisionError,TypeError) as e:
        print(type(e),e)
        return None

if __name__=='__main__':
    print(parseField0PolyFromStr("x**4+4*x**5+2*x+2"))
    
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,1]),"T_1") # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    print(parseExpressionFromStr("x**2/(x+1) * T_1", FE.fieldTower))