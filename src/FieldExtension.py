'''
Created on 22.09.2018

@author: Leonard
'''

#Field extension types

TRANS_LOG = 0
TRANS_EXP = 1
ALGEBRAIC = 2

BASE_FIELD = 0 # C
##BASEFUNCTION_FIELD = 1 # field of rational functions in x: C(x)
#EXTENDED_FIELD = 2 # extended field: C(x,T)
BASE_VARIABLE = "x"
EXTENSION_VARIABLE = "T"
#VARIABLES = [BASE_VARIABLE]#,EXTENSION_VARIABLE]

#fieldExtension = None

fieldTower = None

"""def updateVariables():
    if fieldTower==None:
        return
    if fieldTower.towerHeight>=len(VARIABLES):
        for i in range(len(VARIABLES),fieldTower.towerHeight+1):
            VARIABLES.append("T_{{{}}}".format(i))
            """
class FieldExtension(object):
    '''
    classdocs
    '''

    def __init__(self, extensionType, characteristicFunction, variable):
        """
        logarithmic extension: f' = u'/u where u=characteristicFuntion
        exponential extension f'=u'f where u=characteristicFuntion
        u is in the field extension, that is one level lower
        variable : name of the extensionVariable, usually T (for Theta which is often used in literature)
        """
        self.extensionType = extensionType
        if (extensionType==ALGEBRAIC):
            raise NotImplementedError("Algebraic field extensions are not implemented (yet).")
        self.characteristicFunction = characteristicFunction # = u
        self.variable = variable
        
    def __str__(self):
        var = "exp" if self.extensionType==TRANS_EXP else "log"
        #out = "C(x,T), T={}({})".format(var,str(self.characteristicFunction)) 
        out = "K({}), {}={}({})".format(self.variable,self.variable,var,str(self.characteristicFunction)) 
        return out
    def strFunc(self):
        var = "exp" if self.extensionType==TRANS_EXP else "log"
        return "{}({})".format(var,str(self.characteristicFunction))
class FieldTower(object):
    def __init__(self, fieldExtension=None,fieldExtensions=None):
        if fieldExtensions!=None:
            self.fieldExtensions = fieldExtensions
        elif fieldExtension!=None:
            self.fieldExtensions = [fieldExtension]
        else:
            self.fieldExtensions = []
            
        #updateVariables()
    
    @property
    def towerHeight(self):
        return len(self.fieldExtensions)
    
    def getFieldExtension(self, index):
        return self.fieldExtensions[index]
    def getLastExtension(self):
        if self.towerHeight==0:
            return None
        return self.getFieldExtension(self.towerHeight-1)
    def getStrippedTower(self, index):
        return FieldTower(fieldExtensions=self.fieldExtensions[0:index])
    
    def addFieldExtension(self, fieldExtension):
        self.fieldExtensions.append(fieldExtension)
        #updateVariables()
        
    def __str__(self):
        out = "C(x,"
        for i in range(self.towerHeight):
            out += self.getFieldExtension(i).variable+","
        out = out.strip(",")
        out += ") where "
        for i in range(self.towerHeight):
            ext = self.getFieldExtension(i)
            var = "exp" if ext.extensionType==TRANS_EXP else "log"
            out += "{} = {}({}); ".format(self.getFieldExtension(i).variable,var,str(ext.characteristicFunction))
        out = out.strip("; ")
        return out