from sympy import *
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

var ('n,z,beta,alpha,X0', positive=True)
#Transfer function
k=beta-I*alpha
hp=Sum(X0*(1-X0)**n*exp(-n*I*k*z),(n, 1, oo)).doit().simplify()
HP=X0*(1-X0)*exp(-I*k*z)/(1-(1-X0)*exp(-I*k*z))
print(hp**2)
print(hp.subs(beta,1).subs(X0,0.1).subs(alpha,1).subs(z,1).evalf())
print(HP.subs(beta,1).subs(X0,0.1).subs(alpha,1).subs(z,1).evalf())