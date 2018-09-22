'''
Created on 22.09.2018

@author: Leonard
'''
from FieldExtension import *
from Polynomial import *
from RationalFunction import RationalFunction

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
    monomials = poly_raw.split("+")
    poly = Polynomial()
    for mon in monomials:
        l = mon.split("**")
        i = complex(0,1)
        power = 0
        if (len(l)>1):
            power = eval(l[1])
        elif l[0].endswith(var):
            power = 1
        coeff_raw = l[0].strip(var).strip("*")
        coeff = 1
        if len(coeff_raw)!=0:
            coeff = eval(coeff_raw)
        poly += Monomial(power,coeff)
    return poly
def fieldExtensionInput():
    f_e_type1 = getInput("Type of field extension (Transcendental: t, Algebraic: a):", ["t","a"])
    if f_e_type1=="a":
        raise NotImplementedError()
    f_e_type2 = getInput("Type of transcendental extension (Exponential: e, Logarithmic: l):",["e","l"])
    f_e_type = TRANS_LOG if f_e_type2=="l" else TRANS_EXP
    bfunc = "exp" if f_e_type == TRANS_EXP else "log"
    print("Field extension will be of the form: C(x,T), where T={}(u) with u element C(x)".format(bfunc))
    print("u(x)=p(x)/q(x)")
    p = inputPolynomial("x", "p")
    q = inputPolynomial("x", "q")
    u = RationalFunction(p,q)
    print("u(x)={}".format(u))
    return FieldExtension(f_e_type,u)
    
def Main():
    fieldExtension = fieldExtensionInput()
    print(fieldExtension)
    while (True):
        print("integrate f, f(T)=p(T)/q(T), p,q element C(x)[T]")
if __name__ == '__main__':
    Main()