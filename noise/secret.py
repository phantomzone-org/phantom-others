from estimator import *

from sage.all import RealField

RR = RealField(256)

class Secret():
    def __init__(self, distr: ND, dimension: int):
        self.distr = distr
        self.dimension = dimension

    def ErrorDistribution(N: int):
        return Secret(distr=ND.DiscreteGaussian(3.19), dimension=N)
    
    def TernarySecret(N: int):
        return Secret(distr=ND.SparseTernary(n=N, p=int(N/4)), dimension=N)
    
    def variance(self):
        return RR(self.distr.stddev) * RR(self.distr.stddev)