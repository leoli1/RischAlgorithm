'''
Created on 23.10.2018

@author: Leonard
'''
import PrettyPrinting as PP

R_UP = u"\N{RIGHT PARENTHESIS UPPER HOOK}"
R_LO = u"\N{RIGHT PARENTHESIS LOWER HOOK}"
L_UP = u"\N{LEFT PARENTHESIS UPPER HOOK}"
L_LO = u"\N{LEFT PARENTHESIS LOWER HOOK}"
R_VERT = u"\N{RIGHT PARENTHESIS EXTENSION}"
L_VERT = u"\N{LEFT PARENTHESIS EXTENSION}"
T_INT = u"\N{TOP HALF INTEGRAL}"
B_INT = u"\N{BOTTOM HALF INTEGRAL}"
VERT_INT = u"\N{INTEGRAL EXTENSION}"
INT = u"\N{INTEGRAL}"
DIAG_DOWN = u"\N{BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT}"
DIAG_UP = u"\N{BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT}"
LOW_LINE = u"\N{LOW LINE}"

FRACTION_MID = u"\N{BOX DRAWINGS LIGHT HORIZONTAL}"

SUMM_TOP = u"\N{SUMMATION TOP}"
SUMM_BOT = u"\N{SUMMATION BOTTOM}"

def bigSigma():
    sm = PP.StringMatrix(data=[[SUMM_TOP],[SUMM_BOT]])
    sm.mainRow = 1
    return sm
def LeftBracket1(height):
    if height==1:
        return ["("]
    t = [L_UP]
    for _ in range(1,height-1):
        t.append(L_VERT)
    return t+[L_LO]
def RightBracket1(height):
    if height==1:
        return [")"]
    t = [R_UP]
    for _ in range(1,height-1):
        t.append(R_VERT)
    return t+[R_LO]

def IntegralSymbol(height):
    if height==1:
        return [INT]
    t = [T_INT]
    for _ in range(1,height-1):
        t.append(VERT_INT)
    return t+[B_INT]

def UpDiag(height):
    sm = PP.StringMatrix(width=height,height=height)
    for i in range(height):
        sm.setChar(i, height-1-i, DIAG_UP)
    return sm

if __name__=="__main__":
    for (l,r,i) in zip(LeftBracket1(4),RightBracket1(4),IntegralSymbol(4)):
        print l,r,i
    bigSigma().pprint()