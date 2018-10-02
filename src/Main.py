'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
from Polynomial import *
from RationalFunction import RationalFunction
from RischAlgorithm import *

from Parse import parseField0PolyFromStr,parseExpressionFromStr

def getInput(msg, outputs):
    e = True
    o=""
    while e:
        o = raw_input(msg).lower()
        if o in outputs:
            e = False
    return o
    
def fieldTowerInput():
    height = input("Number of field extensions? ")
        
    fieldTower = FE.FieldTower()
    for i in range(height):
        print("Field extension of {}".format(fieldTower))
        f_e_type1 = getInput("Type of field extension (Transcendental: t, Algebraic: a):", ["t","a"])
        if f_e_type1=="a":
            raise NotImplementedError()
        
        f_e_type2 = getInput("Type of transcendental extension (Exponential: e, Logarithmic: l):",["e","l"])
        if f_e_type2=="e":
            raise NotImplementedError()
        
        f_e_type = FE.TRANS_LOG if f_e_type2=="l" else FE.TRANS_EXP
        bfunc = "exp" if f_e_type == FE.TRANS_EXP else "log"
        print("Field extension: {}(u)".format(bfunc))
        argFunc_str = raw_input("u = ")
        argFunc = parseExpressionFromStr(argFunc_str, fieldTower)
        fe = FE.FieldExtension(f_e_type,argFunc,"T_{}".format(i+1))
        fieldTower.addFieldExtension(fe)
        
    return fieldTower
    
def Main():
    FE.fieldTower = fieldTowerInput()
    print(FE.fieldTower)
    if FE.fieldTower.towerHeight>0:
        print("Use field-extension Variables, instead of functions, e.g. f=x*{} instead of f=x*{}".format(FE.fieldTower.getFieldExtension(0).variable, FE.fieldTower.getFieldExtension(0).strFunc()))
    while (True):
        f_str = raw_input("integrate f = ")
        f = parseExpressionFromStr(f_str, FE.fieldTower)
        print("Integral({}) = {}".format(f.printFull(),printIntegral(f, FE.fieldTower)))

        
        
if __name__ == '__main__':
    Main()