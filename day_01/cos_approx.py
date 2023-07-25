#!/usr/bin/env python
"""Space 477: Python: I

cosine approximation function
"""
__author__ = 'Wyatt Spies'
__email__ = 'wyatt.spies@colorado.edu'

from math import factorial
from math import pi


def cos_approx(x, accuracy=10):
    """
    This is a function designed to approximate the cosine function
    with a taylor expansion
    """

    cos_n = [(-1)**n * x**(2*n)/factorial(2*n) for n in range(accuracy+1)]      #gives array of n terms
    cos = sum(cos_n)                                                            #sums array
    return cos



# Will only run if this is run from command line as opposed to imported
if __name__ == '__main__':  # main code block
    print("cos(0) = ", cos_approx(0))
    print("cos(pi) = ", cos_approx(pi))
    print("cos(2*pi) = ", cos_approx(2*pi))
    print("more accurate cos(2*pi) = ", cos_approx(2*pi, accuracy=50))
    assert cos_approx(0) < 1+1.e-2 and cos_approx(0) > 1-1.e-2, "cos(0) is not 1"

