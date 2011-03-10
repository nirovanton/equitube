import matplotlib.pyplot as plt
import numpy as np
import math
import random

# I need to create a class to
def create_line():

    # The random stuff
    m = random.randint(-20,20)
    l = random.uniform(5,12)
    x_cm, y_cm = random.uniform(10,15)
    
    # May need to work a bit on theta for the van der waals approx.
    theta = np.arctan(math.sqrt(m*m))
    x_min = x_cm - l*np.cos(theta)
    x_max = x_cm + l*np.cos(theta)
    y_min = y_cm - l*np.sin(theta)
    y_max = y_cm + l*np.sin(theta)
    b = y_cm - m*x_cm 

    # y_cm = m*x_cm + b 
    # by setting the line equation for each node equal to one another
    # we can solve for the intersecting quardinants. let y = arbitrary,
    # and check that the corresponding x_intercept coordinate is in
    # the range x_min <= x_int <= x_max  if yes: intersect

    return m, l, theta, x_min, x_cm, x_max, y_min, y_cm, y_max

def intersect(x_min, x_max, b1, b2, m1, m2)
    
    # This function evaluates 2 node to see if they intersect within
    # the node range x_min and x_max both need to come from the same
    # node, it doesn't which node.
 
    x_int = (b2-b1)/(m1-m2)
    if x_min <= x_int and x_int <= x_max:
        return True
    return False



'''
fig = plt.figure()
ax = fig.add_subplot(111)
x, y = np.random.rand(2, 2)
x1,y1 = np.random.rand(2, 2)
line = ax.plot(x, y, 'bs-')
line1 = ax.plot(x1, y1, 'bs-', picker=5)

plt.draw()
plt.show()
'''
