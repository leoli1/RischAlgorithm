'''
Created on 28.09.2018

@author: Leonard
'''
from Utils import isNumber,objEqualsNumber
import FieldExtension as FE
import Polynomial as Pol
import FieldTowerStructure as FTS
import Number

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
        
    def getLogExpression(self, argFunction):
        for log in self.logExpressions:
            if log.argFunction == argFunction:
                return log
        return None
    def getNewLogExpressionsInFieldTower(self, fieldTower,fullTower):
        """
        returns logarithms appearing in the integral (self) that are not in fieldTower
        e.g. fieldTower = C(x,log(x)) = F
             integral = log(x) + log(log(x))
             will return [log(log(x))] since log(x) el F but log(log(x)) not el F
        """
        logs = []
        self.combine_logs()
        for log in self.logExpressions:
            #if log.argFunction.getFieldTower()==fieldTower: # log must not be in fieldTower, e.g. log(x+1) is not el C(x,log(x)), but log(x^2) is el C(x,log(x)), since log(x^2)=2*log(x)
            if not FTS.logarithmIsInFieldTower(log, fieldTower):
                logs.append(log)
        return logs
    
    def simplify(self):
        self.combine_rationals()
        self.combine_logs()
        
    def combine_rationals(self):
        s = 0
        for r in self.poly_rational_partExpressions:
            s += r
            
        if id(s)!=id(0):
            s = s.simplified()
            
        if len(self.poly_rational_partExpressions)!=0:
            if s==0:
                self.poly_rational_partExpressions = []
            else:
                self.poly_rational_partExpressions = [s]
                
    def combine_logs(self):
        newLogs = []
        for l in self.logExpressions:
            log = next((x for x in newLogs if x.argFunction==l.argFunction),None)
            if log==None:
                newLogs.append(l)
            else:
                log.factor += l.factor
        finalLogs = []
        for l in newLogs:
            if l.factor!=0:
                finalLogs.append(l)
        self.logExpressions = finalLogs
        
    def asFunction(self):
        func = 0
        for p in self.poly_rational_partExpressions:
            func += p
        #if self.logExpressions!=[]:
        #    raise NotImplementedError()
        newFieldTower = FE.fieldTower.copy()
        for log in self.logExpressions:
            tower = FE.hasFieldExtension(FE.TRANS_LOG, log.argFunction, FE.fieldTower)
            if tower==None:
                newFieldTower = newFieldTower.addFieldExtension(FE.FieldExtension(FE.TRANS_LOG,log.argFunction,FE.Variable("a"+str(FE.tempVariableNum))))
                newTower = newFieldTower.copy()
                logExpr = Pol.Polynomial([0,log.factor], newFieldTower.getLastVariable())
                func += logExpr
                #raise NotImplementedError("")
            else:
                logExpr = Pol.Polynomial([0,log.factor], tower.getLastVariable())
                func += logExpr
        
        if self.rootSums!=[]:
            raise NotImplementedError()
        
        return func
    def asFunctionWithRootSums(self):
        rs = self.rootSums
        self.rootSums=[]
        f = self.asFunction()
        self.rootSums=rs
        return [f]+rs
    
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
    def __repr__(self):
        return self.__str__()
    
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
    
class NonElementaryIntegral(Integral):
    def __str__(self):
        return "Integral({})".format(super(NonElementaryIntegral,self).__str__())
    def __repr__(self):
        return self.__str__()
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
                if type(self.factor)==Number.SqrRootPolynomial and self.factor.b!=0:
                    factor = "({})*".format(self.factor)
                elif self.factor>=0:
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