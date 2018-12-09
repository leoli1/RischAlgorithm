'''
Created on 22.09.2018

@author: Leonard
'''

#Field extension types

TRANS_LOG = 0
TRANS_EXP = 1
ALGEBRAIC = 2
TRANSCENDENTAL_SYMBOL = 3

fieldTower = None

if not "variables" in dir():
    variables = []

if not "tempVariableNum" in dir():
    tempVariableNum = 0

def hasFieldExtension(type,u,tower):
    for i in range(tower.towerHeight):
        ext = tower.getFieldExtension(i)
        if type==ext.extensionType and u==ext.argFunction:
            return tower.getStrippedTower(i+1)
    return None


class Variable(object):
    def __init__(self, stringRepr):
        self.stringRepr = stringRepr
        self.fieldExtension = None
        variables.append(self)
        global tempVariableNum
        tempVariableNum += 1
        
    def __eq__(self, other):
        if type(other)!=Variable:
            return False
        return self.stringRepr==other.stringRepr
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __str__(self):
        return self.stringRepr
    def __repr__(self):
        return self.__str__()
    
    
def getVariable(stringReprs):
    for var in variables:
        if var.stringRepr == stringReprs:
            return var
def Variables(l):
    vs = [] 
    for v in l:
        vs.append(Variable(v))
    return vs

    
    
class FieldExtension(object):
    '''
    stores the data of a particular differential field extension
    '''

    def __init__(self, extensionType, argFunction, variable):
        """
        logarithmic extension: f' = u'/u where u=argFunction
        exponential extension f'=u'f where u=argFunction
        u is in the field extension, that is one level lower
        variable : instance of Variable to represent the field extension
        """
        self.extensionType = extensionType
        if (extensionType==ALGEBRAIC):
            raise NotImplementedError("Algebraic field extensions are not implemented (yet).")
        self.argFunction = argFunction # = u
        
        if type(variable)!=Variable:
            raise TypeError()
        self.variable = variable
        self.variable.fieldExtension = self
        self.fieldTower = None
        
        
    def __eq__(self, other):
        if type(other)!=FieldExtension:
            return False
        if id(self)==id(other):
            return True
        return self.extensionType==other.extensionType and self.argFunction==other.argFunction
    def __ne__(self, other):
        return not self.__eq__(other)
    def getFVar(self):
        return "exp" if self.extensionType==TRANS_EXP else ("log" if self.extensionType==TRANS_LOG else "")
    def __str__(self):
        var = self.getFVar()
        out = "K({}), {}={}({})".format(self.variable,self.variable,var,str(self.argFunction)) 
        return out
    def __repr__(self):
        return self.__str__()
    def strFunc(self):
        var = "exp" if self.extensionType==TRANS_EXP else "log"
        return "{}({})".format(var,self.argFunction.printFull())
    
class FieldTower(object):
    """
    Tower of field extensions
    """
    def __init__(self, fieldExtension=None,fieldExtensions=None,_withoutBaseExtension=False):
        """
        FieldTower() = C(x)
        FieldTower(fieldExtension) = C(x,fieldExtension.variable)
        etc...
        """
        k = []
        if not _withoutBaseExtension:
            k = [BASEEXTENSION]
        if fieldExtensions!=None:
            self.fieldExtensions = k + fieldExtensions
            
        elif fieldExtension!=None:
            self.fieldExtensions = k+[fieldExtension]
        else:
            self.fieldExtensions = k
            
        for fe in self.fieldExtensions:
            if fe.fieldTower==None:
                fe.fieldTower = self.copy()
            

    
    @property
    def towerHeight(self):
        return len(self.fieldExtensions)
    
    def getFieldExtension(self, index):
        return self.fieldExtensions[index]
    def getLastExtension(self):
        if self.towerHeight==0:
            return None
        return self.getFieldExtension(self.towerHeight-1)
    def getLastVariable(self):
        if self.towerHeight==0:
            return BASEVARIABLE
        return self.getLastExtension().variable
    
    def getStrippedTower(self, index):
        ft = FieldTower(_withoutBaseExtension=True)
        ft.fieldExtensions = self.fieldExtensions[0:index]
        return ft
    
    def prevTower(self):
        return self.getStrippedTower(self.towerHeight-1)
    
    def addFieldExtension(self, fieldExtension):
        ft = self.copy()
        ft.fieldExtensions.append(fieldExtension)
        if fieldExtension.fieldTower==None:
            fieldExtension.fieldTower = ft
        return ft
    
    def getLogarithms(self):
        logs = []
        for fext in self.fieldExtensions:
            if fext.extensionType==TRANS_LOG:
                logs.append(fext)
                
        return logs
    
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
        if id(self)==id(other):
            return True
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
        out = "C("
        for i in range(self.towerHeight):
            out += self.getFieldExtension(i).variable.stringRepr+","
        out = out.strip(",")
        out += ")"
        if self.towerHeight==0:
            raise Exception()
        if self.towerHeight==1:
            return out
        out += " where "
        for i in range(self.towerHeight):
            ext = self.getFieldExtension(i)
            if ext.argFunction==None:
                continue
            var = ext.getFVar()
            if i==1:
                out += "{} = {}({}); ".format(ext.variable.stringRepr,var,str(ext.argFunction))
            else:
                out += "{} = {}({}) = {}({}); ".format(ext.variable.stringRepr,var,str(ext.argFunction),var,ext.argFunction.printFull())
        out = out.strip("; ")
        return out.strip(" where ")
    def __repr__(self):
        return str(self)
    
    
if "BASEVARIABLE" not in dir():
    BASEVARIABLE = Variable('x') # base variable for C(x) functions
    BASEEXTENSION = FieldExtension(TRANSCENDENTAL_SYMBOL,None,BASEVARIABLE)
    BASEFIELD = FieldTower(BASEEXTENSION,_withoutBaseExtension=True)
    