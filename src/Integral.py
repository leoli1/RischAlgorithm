'''
Created on 28.09.2018

@author: Leonard
'''
from Utils import isNumber,objEqualsNumber
import FieldExtension as FE
import Polynomial as Pol

class Integral(object):

    def __init__(self,poly_rationals=None,logs=None, rootSums=None):
        if (poly_rationals!=None and type(poly_rationals)!=list) or (logs!=None and type(logs)!=list) or (rootSums!=None and type(rootSums)!=list):
            raise Exception("parameters have to be either None or lists")
        self.poly_rational_partExpressions = [] if poly_rationals == None else poly_rationals
        self.logExpressions = [] if logs==None else logs
        self.rootSums = [] if rootSums==None else rootSums
        
        newPRs = []
        for pr in self.poly_rational_partExpressions:
            if not pr==0 and not pr.isZero():
                newPRs.append(pr)
        self.poly_rational_partExpressions = newPRs
        
    def getNewLogExpressionsInFieldTower(self, fieldTower,fullTower):
        logs = []
        for log in self.logExpressions:
            #if log.argFunction.getFieldTower()==fieldTower: # log must not be in fieldTower, e.g. log(x+1) is not el C(x,log(x)), but log(x^2) is el C(x,log(x)), since log(x^2)=2*log(x)
            if log.argFunction == fullTower.getLastExtension().characteristicFunction:
                logs.append(log)
            else:
                if fieldTower.towerHeight>1:
                    raise NotImplementedError()
                elif fieldTower.towerHeight==0:
                    logs.append(log)
                else:
                    #if not (fieldTower.getLastExtension().characteristicFunction/log.argFunction).isConstant():
                    #    logs.append(log)
                    (q,r) = Pol.PolyDiv(log.argFunction, fieldTower.getLastExtension().characteristicFunction)
                    if r==0 or r.isZero():
                        pass
                    else:
                        logs.append(log)
        return logs
    def combine_rationals(self):
        s = 0
        for r in self.poly_rational_partExpressions:
            s += r
        self.poly_rational_partExpressions = [s]
    def asFunction(self):
        func = 0
        for p in self.poly_rational_partExpressions:
            func += p
        #if self.logExpressions!=[]:
        #    raise NotImplementedError()
        for log in self.logExpressions:
            tower = FE.hasFieldExtension(FE.TRANS_LOG, log.argFunction, FE.fieldTower)
            if tower==None:
                raise NotImplementedError("")
            else:
                logExpr = Pol.Polynomial([0,log.factor], tower)
                func += logExpr
        
        if self.rootSums!=[]:
            raise NotImplementedError()
        
        return func
    
    def __add__(self, other):
        newLogs = self.logExpressions + other.logExpressions
        newPolys = self.poly_rational_partExpressions+other.poly_rational_partExpressions
        rootSums = self.rootSums + other.rootSums
        return Integral(newPolys,newLogs,rootSums)
    
    def __str__(self):
        out = ""
        for a in self.poly_rational_partExpressions:
            out += str(a)+"+"
        for a in self.logExpressions:
            out += str(a)+"+"
        for a in self.rootSums:
            out += str(a)+"+"
        return out.strip("+")
    
    def printFull(self): # replaces fieldextension variables with their functions, i.e. T = log(x)
        out = ""
        for a in self.poly_rational_partExpressions:
            if isNumber(a):
                out += str(a) + "+"
            else:
                out += a.printFull()+"+"
        for a in self.logExpressions:
            out += a.printFull()+"+"
        for a in self.rootSums:
            out += a.printFull()+"+"
        return out.strip("+")
    

class IntegralNotElementaryError(Exception):
    pass


class LogFunction(object):

    def __init__(self, argFunction, factor):

        self.argFunction = argFunction
        self.factor = factor
    def __str__(self):
        factor = ""
        if objEqualsNumber(self.factor, 0):
            return ""
        if isNumber(self.factor):
            if not objEqualsNumber(self.factor, 1):
                if self.factor>=0:
                    factor = str(self.factor)+"*"
                else:
                    factor = "({})*".format(self.factor)
        elif self.factor.isConstant():
            c = self.factor.getConstant()
            if not objEqualsNumber(c, 1):
                if c>=0:
                    factor = str(self.factor.getConstant())+"*"
                else:
                    factor = "({})".format(c)
        else:
            factor = "[{}]*".format(str(self.factor))
        return "{}log({})".format(factor,str(self.argFunction))
    def printFull(self):
        factor = ""
        if isNumber(self.factor):
            if not objEqualsNumber(self.factor, 1):
                if self.factor>=0:
                    factor = str(self.factor)+"*"
                else:
                    factor = "({})*".format(self.factor)
        elif self.factor.isConstant():
            c = self.factor.getConstant()
            if not objEqualsNumber(c, 1):
                if type(c)==complex or c>=0:
                    factor = str(self.factor.getConstant())+"*"
                else:
                    factor = "({})".format(c)
        else:
            factor = "[{}]".format(self.factor.printFull())
        return "{}log({})".format(factor,self.argFunction.printFull())