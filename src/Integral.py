'''
Created on 28.09.2018

@author: Leonard
'''

class Integral(object):

    def __init__(self, content):
        self.content = str(content)
        self.poly_rational_partExpressions = None
        self.logExpressions = None
        
    def __add__(self, other):
        if len(other.content)==0:
            return self
        elif len(self.content)==0:
            return other
        return Integral("{}+{}".format(str(self), str(other)))
    def __sub__(self, other):
        return Integral("{}-{}".format(str(self), str(other)))
    def __str__(self):
        return self.content
    
    def printFull(self): # replaces fieldextension variables with their functions, i.e. T = log(x)
        out = str(self)
        return out # todo
    

class IntegralNotElementaryException(Exception):
    pass