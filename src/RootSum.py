'''
Created on 28.09.2018

@author: Leonard
'''

class RootSum(object):


    def __init__(self, poly, expr, exprVar="x"):
        self.poly = poly
        self.expr = expr
        self.exprVar = exprVar
        
    def __str__(self):
        if type(self.expr)!=str:
            return "RootSum(w | {}=0, {})".format(self.poly.strCustomVar("w"),self.expr)
        return "RootSum(w | {}=0, {})".format(self.poly.strCustomVar("w"),self.expr.replace(self.exprVar,"w"))
    def printFull(self):
        #TODO
        return str(self)