'''
Created on 22.09.2018

@author: Leonard
'''

#Field extension types

TRANS_LOG = 0
TRANS_EXP = 1
ALGEBRAIC = 2
TRANSCENDENTAL_SYMBOL = 3

BASE_FIELD = 0 # C
##BASEFUNCTION_FIELD = 1 # field of rational functions in x: C(x)
#EXTENDED_FIELD = 2 # extended field: C(x,T)
BASE_VARIABLE = "x"
EXTENSION_VARIABLE = "T"
#VARIABLES = [BASE_VARIABLE]#,EXTENSION_VARIABLE]

#fieldExtension = None

fieldTower = None

def hasFieldExtension(type,u,tower):
    for i in range(tower.towerHeight):
        ext = tower.getFieldExtension(i)
        if type==ext.extensionType and u==ext.characteristicFunction:
            return tower.getStrippedTower(i+1)
    return None


class FieldExtension(object):
    '''
    stores the data of a differential field extension
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
        
        
    def __eq__(self, other):
        if type(other)!=FieldExtension:
            return False
        return self.extensionType==other.extensionType and self.characteristicFunction==other.characteristicFunction
    def __ne__(self, other):
        return not self.__eq__(other)
    def getFVar(self):
        return "exp" if self.extensionType==TRANS_EXP else ("log" if self.extensionType==TRANS_LOG else "")
    def __str__(self):
        var = self.getFVar()
        out = "K({}), {}={}({})".format(self.variable,self.variable,var,str(self.characteristicFunction)) 
        return out
    def __repr__(self):
        return self.__str__()
    def strFunc(self):
        var = "exp" if self.extensionType==TRANS_EXP else "log"
        return "{}({})".format(var,self.characteristicFunction.printFull())
        #return "{}({})".format(var,str(self.characteristicFunction))
    
class FieldTower(object):
    """
    Tower of field extensions
    """
    def __init__(self, fieldExtension=None,fieldExtensions=None):
        """
        FieldTower() = C(x)
        FieldTower(fieldExtension) = C(x,fieldExtension.variable)
        etc...
        """
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
    def prevTower(self):
        return self.getStrippedTower(self.towerHeight-1)
    
    def addFieldExtension(self, fieldExtension):
        self.fieldExtensions.append(fieldExtension)

    def copy(self):
        return self.getStrippedTower(self.towerHeight)
    
    def isExtendedTowerOf(self, other):
        """
        tests if other.fieldExtensions is a subset of self.fieldExtensions
        """
        if other.towerHeight>=self.towerHeight:
            return False
        for i in range(other.towerHeight):
            if self.getFieldExtension(i) != other.getFieldExtension(i):
                return False
        return True
    
    def __eq__(self, other):
        if type(other)!=FieldTower:
            return False
        if self.towerHeight != other.towerHeight:
            return False
        for i in range(self.towerHeight):
            if self.getFieldExtension(i)!=other.getFieldExtension(i):
                return False
            
        return True
    def __ne__(self, other):
        return not (self==other)
    def __str__(self):
        out = "C(x,"
        for i in range(self.towerHeight):
            out += self.getFieldExtension(i).variable+","
        out = out.strip(",")
        out += ")"
        if self.towerHeight==0:
            return out
        out += " where "
        for i in range(self.towerHeight):
            ext = self.getFieldExtension(i)
            var = ext.getFVar()
            if i==0:
                out += "{} = {}({}); ".format(ext.variable,var,str(ext.characteristicFunction))
            else:
                out += "{} = {}({}) = {}({}); ".format(ext.variable,var,str(ext.characteristicFunction),var,ext.characteristicFunction.printFull())
        out = out.strip("; ")
        return out
    def __repr__(self):
        return str(self)
    