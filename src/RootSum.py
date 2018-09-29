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
        
        return "RootSum(w | {}=0, {})".format(self.poly.strCustomVar("w"),self.expr.replace(self.exprVar,"w"))
        