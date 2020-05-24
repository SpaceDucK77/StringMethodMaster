import numpy as np
from md_tools2 import*

tested = [get_extremes, normalise]
altered = [delta_angles]

a = np.array([list(range(11))])
b = np.array(a)
np.random.shuffle(b.T)
c = np.array(a+355)%360-180
d = np.array(b+355)%360-180
A = np.concatenate((a,b,c,d))
Ao = [True]*4

A=A.T
print(A)
mins, deltas = get_extremes(A,Ao)
#print (mins)
print(deltas)
NA = normalise(A, mins, deltas, 4)
#print(NA)
unNA, undr = denormalise(NA, A, mins, deltas, 4)
print(unNA)
print(undr)
