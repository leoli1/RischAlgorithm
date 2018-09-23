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

fieldExtension = None

class FieldExtension(object):
    '''
    classdocs
    '''

    def __init__(self, extensionType, characteristicFunction):
        '''
        logarithmic extension: f' = u'/u where u=characteristicFuntion
        exponential extension f'=u'f where u=characteristicFuntion
        u is in the field extension, that is one level lower
        '''
        self.extensionType = extensionType
        if (extensionType==ALGEBRAIC):
            raise NotImplementedError("Algebraic field extensions are not implemented (yet).")
        self.characteristicFunction = characteristicFunction # = u
        
    def __str__(self):
        var = "exp" if self.extensionType==TRANS_EXP else "log"
        #out = "C(x,T), T={}({})".format(var,str(self.characteristicFunction)) 
        out = "K(T), T={}({})".format(var,str(self.characteristicFunction)) 
        return out
        
class FieldTower(object):
    def __init__(self, fieldExtension=None,fieldExtensions=None):
        if fieldExtensions!=None:
            self.fieldExtensions = fieldExtensions
        elif fieldExtension!=None:
            self.fieldExtensions = [fieldExtension]
        else:
            self.fieldExtensions = []
        
    
    @property
    def towerHeight(self):
        return len(self.fieldExtensions)
    
    def getFieldExtension(self, index):
        return self.fieldExtensions[index]
    def getStrippedExtensions(self, index):
        return FieldTower(fieldExtensions=self.fieldExtensions[0:index])
    
    def addFieldExtension(self, fieldExtension):
        self.fieldExtensions.append(fieldExtension)
        
    def __str__(self):
        out = "C("
        for i in range(self.towerHeight+1):
            out += VARIABLES[i]
        out = out.strip(",")
        out += ") where"
        for i in range(self.towerHeight):
            ext = self.getFieldExtension(i)
            var = "exp" if ext.extensionType==TRANS_EXP else "log"
            out += "{} = {}({}); ".format(VARIABLES[i+1],var,str(ext.characteristicFuntion))
        out = out.strip(";")
        return out