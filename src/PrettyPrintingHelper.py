'''
Created on 23.10.2018

@author: Leonard
'''

R_UP = u"\N{RIGHT PARENTHESIS UPPER HOOK}"
R_LO = u"\N{RIGHT PARENTHESIS LOWER HOOK}"
L_UP = u"\N{LEFT PARENTHESIS UPPER HOOK}"
L_LO = u"\N{LEFT PARENTHESIS LOWER HOOK}"
R_VERT = u"\N{RIGHT PARENTHESIS EXTENSION}"
L_VERT = u"\N{LEFT PARENTHESIS EXTENSION}"

FRACTION_MID = u"\N{BOX DRAWINGS LIGHT HORIZONTAL}"

def LeftBracket1(height):
    if height==1:
        return ["("]
    t = [L_UP]
    for i in range(1,height-1):
        t.append(L_VERT)
    return t+[L_LO]
def RightBracket1(height):
    if height==1:
        return [")"]
    t = [R_UP]
    for i in range(1,height-1):
        t.append(R_VERT)
    return t+[R_LO]


if __name__=="__main__":
    for c in (LeftBracket1(4)):
        print c