'''
Created on 22.09.2018

@author: Leonard
'''
import FieldExtension as FE
from Polynomial import *
from RationalFunction import RationalFunction
from RischAlgorithm import *

from Parse import parseField0PolyFromStr

def getInput(msg, outputs):
    e = True
    o=""
    while e:
        o = raw_input(msg).lower()
        if o in outputs:
            e = False
    return o
def inputPolynomial(var, pname):
    poly_raw = raw_input("{}({})=".format(pname,var)) # polynomial = sum of the monomials separated by '+'
    return parseField0PolyFromStr(poly_raw,var=var)

def inputRational(var,pnames):
    p = inputPolynomial("x", pnames[0])
    q = inputPolynomial("x", pnames[1])
    return RationalFunction(p,q)
def inputPolynomialFromExtendedField(var, pname):
    deg = input("deg({})=".format(pname))
    poly_raw = ""
    for i in range(deg,-1,-1):
        poly_raw += "a_{{{}}}{}**{}+".format(i,var,i)
    poly_raw = poly_raw.strip("+")
    print(poly_raw)
    print("a_i element C[x]") # C[x] instead of C(x) since you can multiply both p(T) and q(T) with a common multiple of the denominators in their coefficients
    coeffs = []
    for i in range(deg,-1,-1):
        coeffs.append(inputPolynomial("x","a_{}".format(i)))
    coeffs.reverse()
    return Polynomial(coefficients=coeffs,field=1)
def fieldExtensionInput(fieldTower=None):
    f_e_type1 = getInput("Type of field extension (Transcendental: t, Algebraic: a):", ["t","a"])
    if f_e_type1=="a":
        raise NotImplementedError()
    f_e_type2 = getInput("Type of transcendental extension (Exponential: e, Logarithmic: l):",["e","l"])
    f_e_type = FE.TRANS_LOG if f_e_type2=="l" else FE.TRANS_EXP
    bfunc = "exp" if f_e_type == FE.TRANS_EXP else "log"
    print("Field extension will be of the form: C(x,T), where T={}(u) with u element C(x)".format(bfunc))
    print("u(x)=p(x)/q(x)")
    u = inputRational("x",("p","q"))
    print("u(x)={}".format(u))
    return FE.FieldExtension(f_e_type,u,"T_{{0}}")
    
"""def fieldTowerInput():
    height = input("Number of field extensions? ")
    FE.VARIABLES = [FE.BASE_VARIABLE]
    for i in range(height):
        FE.VARIABLES.append("T_{{{}}}".format(i+1))
        
    fieldTower = FE.FieldTower()
    for i in range(height):
        t = "C(x,"
        for j in range(1,i+2):
            t+= FE.VARIABLES[j]+","
        t = t.strip(",") + ")"
        print(t)
        
    return fieldTower"""
    
def Main():
    fieldExtension = fieldExtensionInput()
    FE.fieldTower = FE.FieldTower(fieldExtension=fieldExtension)
   # FE.updateVariables()
    print(FE.fieldTower)
    while (True):
        print("integrate f, f(T)=p(T)/q(T); p,q element C(x)[T]")
        print("input p(T)")
        p = inputPolynomialFromExtendedField("T", "p")
        print("input q(T)")
        q = inputPolynomialFromExtendedField("T", "q")
        f = RationalFunction(p,q,field=1)
        print(f)
        #print(f.differentiate())
        output = Integrate(f, FE.fieldTower)
        if output!=None:
            print("Integral({} dx) = {} + C".format(str(f),output))
        else:
            print("Integral({} dx) is not elementary".format(str(f)))
        
        
if __name__ == '__main__':
    Main()