# skew_norm.py
#
# Functions to define a skew normal distribution from
# SciPy.  1, 2, and 3 component are defined

from scipy.stats import norm

def pdf(x, x0=0.0, sigma=1.0, alpha=0.0):
    arg = (x - x0) / sigma
    return  2.0 / sigma * norm.pdf(arg) * norm.cdf(alpha * arg)

def pdf_2comp(x, x1, s1, a1, x2, s2, a2, w):
    return pdf(x, x1, s1, a1) + (1.0 - w) * pdf(x, x2, s2, a2)

def pdf_3comp(x, x1, s1, a1, w1, x2, s2, a2, w2, x3, s3, a3, w3):
    return 1.0 / (w1 + w2 + w3) * (w1 * pdf(x, x1, s1, a1)
                                 + w2 * pdf(x, x2, s2, a2)
                                 + w3 * pdf(x, x3, s3, a3))


