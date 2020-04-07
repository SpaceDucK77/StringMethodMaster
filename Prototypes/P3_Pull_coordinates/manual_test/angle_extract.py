from md_tools import *
import time

def get_angles(path = ""):
    angles = []
    for i in range(501):
        angles.append(get_angle(path+"conf" + str(i) + ".gro", "test.ndx"))
        #time.sleep(.01)
    return(angles)

for angle in get_angles():
    print(angle)
