'''
Created on 23.09.2018

@author: Leonard
'''

class LogFunction(object):
    '''
    classdocs
    '''


    def __init__(self, argFunction):
        '''
        Constructor
        '''
        self.argFunction = argFunction
        self.factor = 1
        self.summand = 0
    def differentiate(self):
        return self.argFunction.differentiate()/self.argFunction
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self,other):
        if other==0:
            return self
        self.summand += other
        return self
    def __rmul__(self, other):
        return self.__mul(other)
    def __mul__(self, other):
        if other==1:
            return self
        self.factor *= other
        return self
    def __str__(self):
        return "[{}]+[{}]*log({})".format(self.summand, self.factor, self.argFunction)