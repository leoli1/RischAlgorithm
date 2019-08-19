'''
Created on 23.10.2018

@author: Leonard
'''
import PrettyPrinting as PP
from Main import useUnicode

R_UP = u"\N{RIGHT PARENTHESIS UPPER HOOK}" if useUnicode else "\\"
R_LO = u"\N{RIGHT PARENTHESIS LOWER HOOK}" if useUnicode else "/"
L_UP = u"\N{LEFT PARENTHESIS UPPER HOOK}" if useUnicode else "/"
L_LO = u"\N{LEFT PARENTHESIS LOWER HOOK}" if useUnicode else "\\"
R_VERT = u"\N{RIGHT PARENTHESIS EXTENSION}" if useUnicode else "|"
L_VERT = u"\N{LEFT PARENTHESIS EXTENSION}" if useUnicode else "|"
T_INT = u"\N{TOP HALF INTEGRAL}" if useUnicode else "/"
B_INT = u"\N{BOTTOM HALF INTEGRAL}" if useUnicode else "/"
VERT_INT = u"\N{INTEGRAL EXTENSION}" if useUnicode else "|"
INT = u"\N{INTEGRAL}"
DIAG_DOWN = u"\N{BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT}" if useUnicode else "\\"
DIAG_UP = u"\N{BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT}" if useUnicode else "/"
LOW_LINE = u"\N{LOW LINE}" if useUnicode else "_"

FRACTION_MID = u"\N{BOX DRAWINGS LIGHT HORIZONTAL}" if useUnicode else "-"

SUMM_TOP = u"\N{SUMMATION TOP}"
SUMM_BOT = u"\N{SUMMATION BOTTOM}"

def bigSigma():
    if useUnicode:
        sm = PP.StringMatrix(data=[[SUMM_TOP],[SUMM_BOT]])
        sm.mainRow = 1
    else:
        sm = PP.StringMatrix(data=[[" ","_"],
                                   ["\\"," "],
                                   ["/"," "],
                                   [" ","-"]])
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
    if height==1 and useUnicode:
        return [INT]
    if not useUnicode and height<=2:
        height=3
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
        print (l,r,i)
    bigSigma().pprint()