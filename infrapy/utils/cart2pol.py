import math
def cart2pol(x, y):
    '''
    Transform Cartesian to polar coordinates

    Inputs:
    x is the x coordinate
    y is the y coordinate

    Outputs:
    th is the angle
    r is the range
    '''
    th = math.atan2(y, x)
    r = math.hypot(x, y)
    return th, r
