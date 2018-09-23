
# coding: utf-8

# In[17]:


import scipy.integrate as sci
import numpy as np
import math
from matplotlib import pyplot

N = 80
T1 = []
T2 = []
T3 = []
A = np.linspace(0.1,4,N)

for i in A:
    dT = lambda b:((2*2**0.5)/(i**4-b**4)**0.5)
    b = np.linspace(0,i,N,endpoint=False)
    dT1vec = np.vectorize(dT)
    dT1_b = dT1vec(b)
    T1.append(sci.trapz(dT1_b, b)) # calculate time period by Trapezoidal integration method
    T2.append(sci.quad(dT,0,i)[0]) # calculate time period by normal python integration
    T3.append(sci.simps(dT1_b, b)) # calculate time period by Simpson's rule
    
    
pyplot.plot(b,T1,'g*-', label ="Time period by Trapezoidal rule")
pyplot.plot(b,T2,'r--', label ="Time period by normal python integration rule")
pyplot.plot(b,T3,'b--', label ="Time period by Simpson's rule")
pyplot.xlabel('Amplitude')
pyplot.ylabel('Period')
pyplot.title("T vs. A by different integration methods")
pyplot.legend()

