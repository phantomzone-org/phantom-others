from __future__ import annotations

from estimator import *
from sage.all import ceil, exp, log, RealField, sqrt, Integer, erf
from enum import Enum

RR = RealField(256)


def format_rr(v: RR):
    return f"{v.numerical_approx().str()} (log2={v.log2().numerical_approx().str()})"

class Decomposer():
    def __init__(self, d_a: Integer, d_b: Integer, logQ: Integer, logB: Integer):
        assert logQ  >= (d_a * logB)
        if d_b is not None:
            assert logQ  >= (d_b * logB)

        self.ignore_bits_a = logQ - (d_a*logB)
        if self.ignore_bits_a < 0:
            self.ignore_bits_a = 0
        if d_b is not None:
            self.ignore_bits_b = logQ - (d_b*logB)
            if self.ignore_bits_b < 0:
                self.ignore_bits_b = 0
        else:
            self.ignore_bits_b = None

        self.d_a = d_a
        self.d_b = d_b
        self.logB = logB
        self.logQ = logQ
    
    @staticmethod
    def double_decomposer(d_a: Integer, d_b: Integer, logQ: Integer, logB: Integer):
        return Decomposer(d_a=d_a, d_b=d_b, logQ=logQ, logB=logB)
    
    @staticmethod
    def single_decomposer(d: Integer, logQ: Integer, logB: Integer):
        return Decomposer(d_a=d, d_b=None, logB=logB, logQ=logQ)

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

class ParameterVariant(Enum):
    INTERACTIVE_MULTIPARTY = 1
    NON_INTERACTIVE_MULTIPARTY = 2

class Parameters():
    def __init__(
        self,
        logQ, 
        logQ_ks,
        logq,
        logN, 
        n,
        w,
        lwe_sk: Secret,
        rlwe_sk: Secret,
        rgsw_by_rgsw_decomposer: Decomposer,
        rlwe_by_rgsw_decomposer: Decomposer,
        auto_decomposer: Decomposer,
        lwe_decomposer: Decomposer,
        non_interactive_uitos_decomposer: Decomposer,
        fresh_noise_std,
        variant: ParameterVariant,
        parties,
    ):
        self.k = parties


        if variant == ParameterVariant.NON_INTERACTIVE_MULTIPARTY:
            assert non_interactive_uitos_decomposer is not None
            assert non_interactive_uitos_decomposer.d_b==None
            assert non_interactive_uitos_decomposer.logQ == logQ
        else:
            assert non_interactive_uitos_decomposer is None

        # There's only single d for LWE decomposer
        assert lwe_decomposer.d_b == None

        # There's only single d for RLWE Auto decomposer
        assert auto_decomposer.d_b == None

        assert rgsw_by_rgsw_decomposer.d_a is not None and rgsw_by_rgsw_decomposer.d_b is not None
        assert rlwe_by_rgsw_decomposer.d_a is not None and rlwe_by_rgsw_decomposer.d_b is not None

        assert rgsw_by_rgsw_decomposer.logQ == logQ
        assert rlwe_by_rgsw_decomposer.logQ == logQ
        assert auto_decomposer.logQ == logQ
        assert lwe_decomposer.logQ == logQ_ks


        self.rgsw_by_rgsw_decomposer = rgsw_by_rgsw_decomposer
        self.rlwe_by_rgsw_decomposer = rlwe_by_rgsw_decomposer
        self.lwe_decomposer = lwe_decomposer
        self.auto_decomposer = auto_decomposer
        self.non_interactive_uitos_decomposer = non_interactive_uitos_decomposer


        self.logQ = logQ
        self.Q = 1<<logQ

        self.logQ_ks = logQ_ks
        self.Q_ks = 1<<logQ_ks

        self.logq = logq
        self.q = 1<<logq

        self.logN = logN
        self.N = 1<<logN

        self.n = n

        self.fresh_noise_std = RR(fresh_noise_std)
        self.var = RR(fresh_noise_std)*RR(fresh_noise_std)

        self.w = w

        self.variant = variant

        assert lwe_sk.dimension == self.n
        assert rlwe_sk.dimension == self.N

        self.lwe_sk = lwe_sk
        self.rlwe_sk = rlwe_sk


    def noise_multi_party(self):
        n = RR(self.n)
        N = RR(self.N)

        k = RR(self.k) # no. of parties
        
        var = self.var # 3.19^2
        # rlwe/lwe sk variance for k parties
        var_sk_rlwe = k*self.rlwe_sk.variance()
        var_sk_lwe = k*self.lwe_sk.variance()

        # Fresh RGSW encryption
        var_fresh = RR(0)
        match self.variant:
            case ParameterVariant.INTERACTIVE_MULTIPARTY:
                var_fresh = (N * self.rlwe_sk.variance() * k * var) + var*((N*var_sk_rlwe)+1)
            case ParameterVariant.NON_INTERACTIVE_MULTIPARTY:
                B_uitos = RR(1<<self.non_interactive_uitos_decomposer.logB)
                var_fresh = (k * var) * (B_uitos*B_uitos)/RR(12)*RR(self.non_interactive_uitos_decomposer.d_a)*N
                var_fresh += (
                        RR(1<<(self.non_interactive_uitos_decomposer.ignore_bits_a*2))/12   
                    *   N
                    *   self.rlwe_sk.variance()
                )
                var_fresh += (k*var*self.rlwe_sk.variance())

        
        # RGSW x RGSW products
        B_rgsw_rgsw = RR(1<<self.rgsw_by_rgsw_decomposer.logB)
        d_a_rgsw_by_rgsw = self.rgsw_by_rgsw_decomposer.d_a
        d_b_rgsw_by_rgsw = self.rgsw_by_rgsw_decomposer.d_b
        var_rgswbyrgsw_a = (d_a_rgsw_by_rgsw * ((B_rgsw_rgsw*B_rgsw_rgsw)/RR(12)) * var_fresh * N)
        tmp = var_rgswbyrgsw_a
        print(f"RGSW x RGSW part A ks noise std: {format_rr(sqrt(tmp))}")
        # Approximation error induced by ignoring some least signifcant bits. 
        # The variance of ignored bits is (2^{ignored_bits})^2
        var_rgswbyrgsw_a += (
                RR(1 << (self.rgsw_by_rgsw_decomposer.ignore_bits_a*2))/12
            *   var_sk_rlwe
            *   N
        )
        print(f"RGSW x RGSW part A inexact noise std: {format_rr(sqrt(var_rgswbyrgsw_a-tmp))}")
        var_rgswbyrgsw_b = (d_b_rgsw_by_rgsw * ((B_rgsw_rgsw*B_rgsw_rgsw)/RR(12) * var_fresh * N))
        var_rgswbyrgsw_b += (
                RR(1 << (self.rgsw_by_rgsw_decomposer.ignore_bits_b*2))/12
        )
        var_brk = k*(var_rgswbyrgsw_a+var_rgswbyrgsw_b)


        # RLWE x RGSW where RGSW has var_brk error variance
        d_a_rlwe_by_rgsw = self.rlwe_by_rgsw_decomposer.d_a
        d_b_rlwe_by_rgsw = self.rlwe_by_rgsw_decomposer.d_b
        B_rlwe_rgsw = RR(1<<self.rlwe_by_rgsw_decomposer.logB)
        var_rlwe_by_rgsw_a = (d_a_rlwe_by_rgsw * ((B_rlwe_rgsw*B_rlwe_rgsw)/12) * (var_brk) * N) 
        tmp = var_rlwe_by_rgsw_a
        print(f"RLWE x RGSW Part A ks noise std: {format_rr(sqrt(tmp))}")
        var_rlwe_by_rgsw_a += (
                RR(1 << (self.rlwe_by_rgsw_decomposer.ignore_bits_a*2))/12
            *   var_sk_rlwe
            *   N
        )
        print(f"RLWE x RGSW Part A inexact noise std: {format_rr(sqrt(var_rlwe_by_rgsw_a-tmp))}")
        var_rlwe_by_rgsw_b = (d_b_rlwe_by_rgsw * ((B_rgsw_rgsw*B_rgsw_rgsw)/12) * (var_brk) * N) 
        var_rlwe_by_rgsw_b += (
                RR(1 << (self.rlwe_by_rgsw_decomposer.ignore_bits_b*2))/12
        )
        var_rlwe_by_rgsw = (var_rlwe_by_rgsw_a+var_rlwe_by_rgsw_b)
        # print("var_rlwe_by_rgsw_a: ", sqrt(var_rlwe_by_rgsw_a))
        # print("var_rlwe_by_rgsw_b: ", sqrt(var_rlwe_by_rgsw_b))


        # Auto
        B_auto = RR(1<<self.auto_decomposer.logB)
        d_auto = self.auto_decomposer.d_a
        var_auto = ((B_auto*B_auto)/12)*(k*var)*N*d_auto
        var_auto += (
            RR(1 << (self.auto_decomposer.ignore_bits_a*2))/12
            *   var_sk_rlwe
            *   N
        )


        # LWE ksk from rlwe secret to lwe secret
        B_lwe = RR(1<<self.lwe_decomposer.logB)
        d_lwe  = self.lwe_decomposer.d_a
        var_ks = (((B_lwe*B_lwe)/12)*(k*var)*d_lwe*N) 
        tmp = var_ks
        print(f"LWE ks noise std: {format_rr(sqrt(tmp))}")
        var_ks += (N*(
                RR(1 << (self.lwe_decomposer.ignore_bits_a*2))/12
            *   var_sk_rlwe
        ))
        print(f"LWE inexact noise std: {format_rr(sqrt(var_ks-tmp))}")


        # var ms1: Q -> Q_ks 
        var_ms1 = ((N*var_sk_rlwe)+1) * RR(1/12)
        # var ms2: Q_ks -> q (odd mod switch. Set variance of rounding error 4/12)
        var_ms2 = ((n*var_sk_lwe)+1) * RR(4/12)

        Q_sq = RR(self.Q*self.Q)
        Q_ks_sq = RR(self.Q_ks*self.Q_ks)
        q_sq = RR(self.q*self.q)
        
        w = RR(self.w)

        # var_acc = (n*var_rlwe_by_rgsw)+(self.q*var_auto)
        worst_case_autos = (((w - 1)/w)*n)+((1/w)*(self.q>>1))
        var_acc = (n*var_rlwe_by_rgsw)+(var_auto*(worst_case_autos))
    
        print(format_rr(sqrt((q_sq*(2*var_acc))/Q_sq)), format_rr(sqrt((q_sq*(var_ms1+var_ks))/Q_ks_sq)), format_rr(sqrt(var_ms1)), format_rr(sqrt(var_ms2)))

        var_zeta_nand = ((q_sq*(2*var_acc))/Q_sq) + ((q_sq*(var_ms1+var_ks))/Q_ks_sq) + var_ms2
        var_zeta_xor = ((q_sq*(4*var_acc))/Q_sq) + ((q_sq*(var_ms1+var_ks))/Q_ks_sq) + var_ms2


        fail_prob_nand = ((RR(self.q)/8)/sqrt(2*var_zeta_nand)).erfc()
        fail_prob_xor = ((RR(self.q)/8)/sqrt(2*var_zeta_xor)).erfc()

        

        print(f'''
            Worst case autos: {format_rr(worst_case_autos)} 
            var_sk_rlwe: {format_rr(var_sk_rlwe)}   
            var: {format_rr(self.var)}
            std_ms1: {format_rr(sqrt(var_ms1))}
            std_ms2: {format_rr(sqrt(var_ms2))}
            std_ks: {format_rr(sqrt(var_ks))}
            std_fresh: {format_rr(sqrt(var_fresh))}
            std_brk: {format_rr(sqrt(var_brk))}
            std_auto: {format_rr(sqrt(var_auto))}
            std_rlwe_by_rgsw: {format_rr(sqrt(var_rlwe_by_rgsw))}
            std_acc: {format_rr(sqrt(var_acc))}
            std_zeta_nand: {format_rr(sqrt(var_zeta_nand))}
            std_zeta_xor: {format_rr(sqrt(var_zeta_xor))}
            failure probability nand: {format_rr(fail_prob_nand)}
            failure probability xor : {format_rr(fail_prob_xor)}
        ''')

        # if fail_prob_nand != D(0):
        print(f'Failure probability nand log 2: {format_rr(fail_prob_nand.log2())}')
        # if fail_prob_nand != D(0):
        print(f'Failure probability xor log 2: {format_rr(fail_prob_xor.log2())}')

    def security(self):
        # LWE
        lwe = LWE.Parameters(n=self.n, q=(1<<self.logQ_ks), Xs=self.lwe_sk.distr, Xe=ND.DiscreteGaussian(3.19),m=self.n)
        lwe_res = LWE.estimate(lwe, red_cost_model = RC.BDGL16)

        print("LWE Security")
        print(lwe)
        print(lwe_res)
        
        print("")

        rlwe = LWE.Parameters(n=self.N, q=(1<<self.logQ), Xs=self.rlwe_sk.distr, Xe=ND.DiscreteGaussian(3.19),m=self.N)
        rlwe_res = LWE.estimate(rlwe, red_cost_model = RC.BDGL16)

        print("RLWE Security")
        print(rlwe)
        print(rlwe_res)



# Interactive 2P; high bandswidth; fast runtime (2ms faster than I_2P_LB_SR but has key size 116Mib whereas I_2P_LB_SR has key size 99.6MiB)
I_2_HB_FR = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=12,
    logN=11, 
    n=520,
    w=10,
    lwe_sk=Secret.ErrorDistribution(N=520),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=6,
        logQ=54,
        d_a=8,
        d_b=7
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13,
    ),
    non_interactive_uitos_decomposer=None,
    fresh_noise_std=3.19,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=2,
)

# Interactive 2P; low bandiwdth; slow runtime (although not that slow)
I_2_LB_SR = Parameters(
    logQ=54, 
    logQ_ks=15,
    logq=11,
    logN=11, 
    n=580,
    w=10,
    lwe_sk=Secret.TernarySecret(N=580),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=7,
        logQ=54,
        d_a=6,
        d_b=5
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=15,
        d=12,
    ),
    non_interactive_uitos_decomposer=None,
    fresh_noise_std=3.19,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=2,
)

I_4 = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=11,
    logN=11, 
    n=620,
    w=10,
    lwe_sk=Secret.TernarySecret(N=620),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=6,
        logQ=54,
        d_a=7,
        d_b=6
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13,
    ),
    non_interactive_uitos_decomposer=None,
    fresh_noise_std=3.19,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=4,
)

I_8 = Parameters(
    logQ=54, 
    logQ_ks=17,
    logq=12,
    logN=11, 
    n=660,
    w=10,
    lwe_sk=Secret.TernarySecret(N=660),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=5,
        logQ=54,
        d_a=9,
        d_b=8
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=17,
        d=14,
    ),
    non_interactive_uitos_decomposer=None,
    fresh_noise_std=3.19,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=8,
)



NI_2 = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=12,
    logN=11, 
    n=520,
    w=10,
    lwe_sk=Secret.ErrorDistribution(N=520),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=4,
        logQ=54,
        d_a=10,
        d_b=9
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13,
    ),
    non_interactive_uitos_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=54,
        d=50
    ),
    fresh_noise_std=3.19,
    variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
    parties=2,
)



NI_4_HB_FR = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=11,
    logN=11, 
    n=620,
    w=10,
    lwe_sk=Secret.TernarySecret(N=620),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=3,
        logQ=54,
        d_a=13,
        d_b=12
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13
    ),
    non_interactive_uitos_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=54,
        d=50
    ),
    fresh_noise_std=3.19,
    variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
    parties=4,
)

NI_4_LB_SR = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=12,
    logN=11, 
    n=620,
    w=10,
    lwe_sk=Secret.TernarySecret(N=620),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=4,
        logQ=54,
        d_a=10,
        d_b=9
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13
    ),
    non_interactive_uitos_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=54,
        d=50
    ),
    fresh_noise_std=3.19,
    variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
    parties=4,
)

NI_8 = Parameters(
    logQ=54, 
    logQ_ks=17,
    logq=12,
    logN=11, 
    n=660,
    w=10,
    lwe_sk=Secret.TernarySecret(N=660),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=2,
        logQ=54,
        d_a=20,
        d_b=18
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=17,
        d=14
    ),
    non_interactive_uitos_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=54,
        d=50
    ),
    fresh_noise_std=3.19,
    variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
    parties=8,
)



# I_2_HB_FR.noise_multi_party()
# I_2_LB_SR.noise_multi_party()
# I_4.noise_multi_party()
# I_8.noise_multi_party()
# NI_2.noise_multi_party()
# NI_2_NEW.noise_multi_party()
# NI_4.noise_multi_party()
# NI_4_HB_FR.noise_multi_party()
# NI_4_LB_SR.noise_multi_party()
# NI_8.noise_multi_party()

# TRIAL.noise_multi_party()



##########
# EXTRAS #
##########


def odd_mod_switch(): 
    lwe_n = 480
    k = RR(2) # No. of parties
    var_sk = Secret.ErrorDistribution(N=lwe_n).variance()

    # ms
    var_ms = ((RR(lwe_n)*(k*var_sk))+1) * RR(4)/RR(12)
    print(f"var ms: {var_ms}")

# odd_mod_switch()

def ksk_noise():
    logB = 12
    d = 2
    logQ = 55
    logN = 11
    var = RR(3.19*3.19)

    B = RR(1<<logB)
    N = RR(1<<logN)

    var_m = var

    # ksk noise
    ksk_var = (B*B)/12 * var * RR(d) * N
    print("std ksk before approximation", sqrt(ksk_var).log2())
    # approximation noise
    ignore_bits = logQ - (d * logB)
    ksk_var += (
            RR(1 << (ignore_bits*2))/12
        *   var_m
        *   N
    )
    
    # B = RR(1<<7)
    # var = ((B*B)/12)*RR(3.19*3.19)*RR(1<<logN)
    print("std ksk", sqrt(ksk_var).log2())

# ksk_noise()

# def kok():
#     lwe = LWE.Parameters(n=509, q=(1<<16), Xs=ND.DiscreteGaussian(3.19), Xe=ND.DiscreteGaussian(3.19),m=500)
#     lwe_res = LWE.estimate(lwe, red_cost_model = RC.BDGL16)
# kok()

# Non_interactive 2 parties with Low communication with failure probability 2^{-48}
NI_2_FP_2_48 = Parameters(
    logQ=54, 
    logQ_ks=15,
    logq=11,
    logN=11, 
    n=480,
    w=10,
    lwe_sk=Secret.ErrorDistribution(N=480),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=5,
        logQ=54,
        d_a=8,
        d_b=7
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=15,
        d=12,
    ),
    non_interactive_uitos_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=54,
        d=50
    ),
    fresh_noise_std=3.19,
    variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
    parties=2,
)

# 8 party Non-interactive 2^{-40} Failure probability
NI_8_FP_2_40 = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=11,
    logN=11, 
    n=520,
    w=10,
    lwe_sk=Secret.ErrorDistribution(N=520),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=2,
        logQ=54,
        d_a=22,
        d_b=21
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13
    ),
    non_interactive_uitos_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=54,
        d=50
    ),
    fresh_noise_std=3.19,
    variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
    parties=8,
)

# Interactive 8P High commuincation, Faster runtime, Failure probability 2^{-42}
I_8_HB_FR = Parameters(
    logQ=54, 
    logQ_ks=16,
    logq=11,
    logN=11, 
    n=520,
    w=10,
    lwe_sk=Secret.ErrorDistribution(N=520),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=4,
        logQ=54,
        d_a=12,
        d_b=11
    ),
    rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
        logB=17,
        logQ=54,
        d_a=1,
        d_b=1
    ),
    auto_decomposer=Decomposer.single_decomposer(
        logB=24,
        logQ=54,
        d=1
    ),
    lwe_decomposer=Decomposer.single_decomposer(
        logB=1,
        logQ=16,
        d=13,
    ),
    non_interactive_uitos_decomposer=None,
    fresh_noise_std=3.19,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=8,
)



# TWO_MP_PARAMS = Parameters(
#     logQ=55, 
#     logQ_ks=16,
#     logq=11,
#     logN=11, 
#     n=520,
#     w=10,
#     lwe_sk=Secret.ErrorDistribution(N=520),
#     rlwe_sk=Secret.TernarySecret(N=1<<11),
#     rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
#         logB=11,
#         logQ=55,
#         d_a=4,
#         d_b=3
#     ),
#     rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
#         logB=11,
#         logQ=55,
#         d_a=2,
#         d_b=1
#     ),
#     auto_decomposer=Decomposer.single_decomposer(
#         logB=11,
#         logQ=55,
#         d=2
#     ),
#     lwe_decomposer=Decomposer.single_decomposer(
#         logB=1,
#         logQ=16,
#         d=13,
#     ),
#     non_interactive_uitos_decomposer=None,
#     fresh_noise_std=3.19,
#     variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
#     parties=8,
# )

# EIGHT_MP_PARAMS = Parameters(
#     logQ=55, 
#     logQ_ks=16,
#     logq=11,
#     logN=11, 
#     n=520,
#     w=10,
#     lwe_sk=Secret.ErrorDistribution(N=520),
#     rlwe_sk=Secret.TernarySecret(N=1<<11),
#     rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
#         logB=11,
#         logQ=55,
#         d_a=4,
#         d_b=3
#     ),
#     rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
#         logB=11,
#         logQ=55,
#         d_a=2,
#         d_b=1
#     ),
#     auto_decomposer=Decomposer.single_decomposer(
#         logB=11,
#         logQ=55,
#         d=2
#     ),
#     lwe_decomposer=Decomposer.single_decomposer(
#         logB=1,
#         logQ=16,
#         d=13,
#     ),
#     non_interactive_uitos_decomposer=None,
#     fresh_noise_std=3.19,
#     variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
#     parties=8,
# )


# NI_2_NEW = Parameters(
#     logQ=54, 
#     logQ_ks=15,
#     logq=11,
#     logN=11, 
#     n=580,
#     w=10,
#     lwe_sk=Secret.TernarySecret(N=580),
#     rlwe_sk=Secret.TernarySecret(N=1<<11),
#     rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
#         logB=4,
#         logQ=54,
#         d_a=10,
#         d_b=9 
#     ),
#     rlwe_by_rgsw_decomposer=Decomposer.double_decomposer(
#         logB=17,
#         logQ=54,
#         d_a=1,
#         d_b=1
#     ),
#     auto_decomposer=Decomposer.single_decomposer(
#         logB=24,
#         logQ=54,
#         d=1
#     ),
#     lwe_decomposer=Decomposer.single_decomposer(
#         logB=1,
#         logQ=15,
#         d=12,
#     ),
#     non_interactive_uitos_decomposer=Decomposer.single_decomposer(
#         logB=1,
#         logQ=54,
#         d=50
#     ),
#     fresh_noise_std=3.19,
#     variant=ParameterVariant.NON_INTERACTIVE_MULTIPARTY,
#     parties=2,
# )