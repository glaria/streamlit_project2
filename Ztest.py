import numpy as np
import math
import scipy.special as scsp
""" Return (z-score, p-value) """

def z2p(z):
    """From z-score return p-value."""
    return 2*(1- (0.5 * (1 + scsp.erf(abs(z) / math.sqrt(2)))))

def zscore(p1, p2, n1, n2): # p1, p2 proportions
	p = (p1*float(n1) + p2*float(n2))/(float(n1) + float(n2))
	numerator = p1 - p2
	denominator = math.sqrt(p*(1-p)*((1/n1)+(1/n2)))
	return numerator/denominator


p1 = 564
p2 = 120
n1 = 6597
n2 = 1990

def main():
    z = zscore(p1/float(n1), p2/float(n2), float(n1), float(n2))
    print(z,z2p(z))
    return (z,z2p(z))

if __name__ == '__main__':
    main()

