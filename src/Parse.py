'''
Created on 28.09.2018

@author: Leonard
'''
import Polynomial as Pol
import RationalFunction as Rat

def parseField0PolyFromStr(poly_raw,var="x"):
    monomials = poly_raw.split("+")
    poly = Pol.Polynomial()
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
        poly += Pol.Monomial(power,coeff)
    return poly
def parseField0RatFromStr(rat_raw,var="x"):
    parts = rat_raw.split("/")
    if len(parts)!=2:
        raise Exception()
    return Rat.RationalFunction(parseField0PolyFromStr(parts[0], var), parseField0PolyFromStr(parts[1], var))

if __name__=='__main__':
    print(parseField0PolyFromStr("x**4+4*x**5+2*x+2"))