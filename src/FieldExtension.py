'''
Created on 22.09.2018

@author: Leonard
'''

#Field extension types

TRANS_LOG = 0
TRANS_EXP = 1
ALGEBRAIC = 2

BASE_FIELD = 0 # C
BASEFUNCTION_FIELD = 1 # field of rational functions in x: C(x)
EXTENDED_FIELD = 2 # extended field: C(x,T)
BASE_VARIABLE = "x"
EXTENSION_VARIABLE = "T"
VARIABLES = [BASE_VARIABLE,EXTENSION_VARIABLE]

class FieldExtension(object):
    '''
    classdocs
    '''

    def __init__(self, extensionType, characteristicFunction):
        '''
        logarithmic extension: f' = u'/u where u=characteristicFuntion
        exponential extension f'=u'f where u=characteristicFuntion
        u element C(x)
        '''
        self.extensionType = extensionType
        if (extensionType==ALGEBRAIC):
            raise NotImplementedError("Algebraic field extensions are not implemented (yet).")
        self.characteristicFuntion = characteristicFunction # = u
        
    def __str__(self):
        var = "exp" if self.extensionType==TRANS_EXP else "log"
        out = "C(x,T), T={}({})".format(var,str(self.characteristicFuntion))
        return out
        
#class FieldTower(object):
#    pass