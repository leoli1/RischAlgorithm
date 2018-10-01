'''
Created on 28.09.2018

@author: Leonard
'''
from Utils import isNumber,objEqualsNumber

class Integral(object):

    def __init__(self,poly_rationals=None,logs=None, rootSums=None):
        if (poly_rationals!=None and type(poly_rationals)!=list) or (logs!=None and type(logs)!=list) or (rootSums!=None and type(rootSums)!=list):
            raise Exception("parameters have to be either None or lists")
        self.poly_rational_partExpressions = [] if poly_rationals == None else poly_rationals
        self.logExpressions = [] if logs==None else logs
        self.rootSums = [] if rootSums==None else rootSums
        
        newPRs = []
        for pr in self.poly_rational_partExpressions:
            if not pr.isZero():
                newPRs.append(pr)
        self.poly_rational_partExpressions = newPRs
        
    def getLogExpressionsInFieldTower(self, fieldTower):
        logs = []
        for log in self.logExpressions:
            if log.argFunction.getFieldTower()==fieldTower:
                logs.append(log)
                
        return logs
    
    def asFunction(self):
        func = 0
        for p in self.poly_rational_partExpressions:
            func += p
        if self.logExpressions!=[]:
            raise NotImplementedError()
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
            out += a.printFull()+"+"
        for a in self.logExpressions:
            out += a.printFull()+"+"
        for a in self.rootSums:
            out += a.printFull()+"+"
        return out.strip("+")
    

class IntegralNotElementaryException(Exception):
    pass


class LogFunction(object):
    '''
    classdocs
    '''


    def __init__(self, argFunction, factor):
        '''
        Constructor
        '''
        self.argFunction = argFunction
        self.factor = factor
    #    self.summand = 0
    #def differentiate(self):
    #    return self.argFunction.differentiate()/self.argFunction
    #def __radd__(self, other):
    #    return self.__add__(other)
    #def __add__(self,other):
    #    if other==0:
    #        return self
    #    self.summand += other
    #    return self
    #def __rmul__(self, other):
    #    return self.__mul(other)
    #def __mul__(self, other):
    #    if other==1:
    #        return self
    #    self.factor *= other
    #    return self
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
                if c>=0:
                    factor = str(self.factor.getConstant())+"*"
                else:
                    factor = "({})".format(c)
        else:
            factor = "[{}]".format(self.factor.printFull())
        return "{}log({})".format(factor,self.argFunction.printFull())
        #return "[{}]+[{}]*log({})".format(self.summand, self.factor, self.argFunction)