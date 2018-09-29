'''
Created on 28.09.2018

@author: Leonard
'''

class Integral(object):

    def __init__(self, content):
        self.content = str(content)
        self.polynomialPart = None
        self.otherPart = None
        
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
        