from decomposer import Decomposer
from secret import Secret
from enum import Enum
from sage.all import RealField, sqrt, log, ceil
from estimator import *

RR = RealField(256)

class ParameterVariant(Enum):
    INTERACTIVE_MULTIPARTY = 1
    NON_INTERACTIVE_MULTIPARTY = 2

def format_rr(v: RR):
    return f"{v.numerical_approx().str()} (log2={v.log2().numerical_approx().str()})"


class Parameters():
    def __init__(
        self,
        msg,
        padding,
        logQ, 
        logQ_ks,
        logq,
        lwe_sk: Secret,
        rlwe_sk: Secret,
        variant: ParameterVariant,
        parties,
        w=10,
        fresh_noise_std=3.19,
        rgsw_by_rgsw_decomposer: Decomposer=None,
        rlwe_by_rgsw_decomposer: Decomposer=None,
        auto_decomposer: Decomposer=None,
        lwe_decomposer:Decomposer=None,
        non_interactive_uitos_decomposer:Decomposer=None,
    ):
        self.k = parties
        self.msg = msg
        self.padding = padding

        if variant == ParameterVariant.NON_INTERACTIVE_MULTIPARTY:
            assert non_interactive_uitos_decomposer is not None
            assert non_interactive_uitos_decomposer.d_b==None
            assert non_interactive_uitos_decomposer.logQ == logQ
        else:
            assert non_interactive_uitos_decomposer is None

        if rlwe_by_rgsw_decomposer == None:
            self.rlwe_by_rgsw_decomposer = Decomposer.double_decomposer(
                logB=1,
                logQ=logQ,
                d_a=logQ,
                d_b=logQ,
            )
        else:
            assert rlwe_by_rgsw_decomposer.d_a is not None and rlwe_by_rgsw_decomposer.d_b is not None
            assert rlwe_by_rgsw_decomposer.logQ == logQ
            self.rlwe_by_rgsw_decomposer = rlwe_by_rgsw_decomposer

        if rgsw_by_rgsw_decomposer == None:
            self.rgsw_by_rgsw_decomposer = Decomposer.double_decomposer(
                logB=1,
                logQ=logQ,
                d_a=logQ,
                d_b=logQ,
            )
        else:
            assert rgsw_by_rgsw_decomposer.d_a is not None and rgsw_by_rgsw_decomposer.d_b is not None
            assert rgsw_by_rgsw_decomposer.logQ == logQ
            self.rgsw_by_rgsw_decomposer = rgsw_by_rgsw_decomposer

        if auto_decomposer == None:
            self.auto_decomposer = Decomposer.single_decomposer(
                logB=1,
                logQ=logQ,
                d=logQ,
            )
        else:
            assert auto_decomposer.logQ == logQ
            assert auto_decomposer.d_b == None
            self.auto_decomposer = auto_decomposer

        if lwe_decomposer == None:
            self.lwe_decomposer = Decomposer.single_decomposer(
                logB=1,
                logQ=logQ_ks,
                d=logQ_ks,
            )
        else:
            # There's only single d for LWE decomposer
            assert lwe_decomposer.d_b == None
            assert lwe_decomposer.logQ == logQ_ks
            self.lwe_decomposer = lwe_decomposer

        self.non_interactive_uitos_decomposer = non_interactive_uitos_decomposer

        self.logQ = logQ
        self.Q = 1<<logQ

        self.logQ_ks = logQ_ks
        self.Q_ks = 1<<logQ_ks

        self.logq = logq
        self.q = 1<<logq

        self.N = rlwe_sk.dimension
        self.n = lwe_sk.dimension

        self.fresh_noise_std = RR(fresh_noise_std)
        self.var = RR(fresh_noise_std)*RR(fresh_noise_std)

        self.w = w

        self.variant = variant

        self.lwe_sk = lwe_sk
        self.rlwe_sk = rlwe_sk

    # TODO
    def key_size(self) -> int:
        return 0

    # TODO
    def complexity(self) -> int:
        return 0

    def optimize(self, log_fail: float) -> float:

        pr_fail = self.pr_fail()

        if pr_fail > log_fail:
            raise Exception(f'initial PR[fail]={pr_fail} > target PR[fail]={log_fail}')

        # TODO: TEST WITH SWITCHED ORDER AND TAKE BEST BUT REQUIRES WAY TO 
        # EVALUATE THE COMPLEXITY OF PARAMETERS
        pr_fail = self.lwe_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.auto_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.rlwe_by_rgsw_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.rgsw_by_rgsw_decomposer.optimize(self.pr_fail, log_fail)
        
        return pr_fail

    def print_decomposers(self):
        print(f'rgsw_by_rgsw_decomposer: logB={self.rgsw_by_rgsw_decomposer.logB:<2d} d_a={self.rgsw_by_rgsw_decomposer.d_a:<2d} d_b={self.rgsw_by_rgsw_decomposer.d_b:<2d}')
        print(f'rlwe_by_rgsw_decomposer: logB={self.rlwe_by_rgsw_decomposer.logB:<2d} d_a={self.rlwe_by_rgsw_decomposer.d_a:<2d} d_b={self.rlwe_by_rgsw_decomposer.d_b:<2d}')
        print(f'auto_decomposer        : logB={self.auto_decomposer.logB:<2d} d_a={self.auto_decomposer.d_a:<2d}')
        print(f'lwe_decomposer         : logB={self.lwe_decomposer.logB:<2d} d_a={self.lwe_decomposer.d_a:<2d}')

    def pr_fail(self):

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
                        RR(1<<(self.non_interactive_uitos_decomposer.ignore_bits_a()*2))/12   
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
        #print(f"RGSW x RGSW part A ks noise std: {format_rr(sqrt(tmp))}")
        # Approximation error induced by ignoring some least signifcant bits. 
        # The variance of ignored bits is (2^{ignored_bits})^2
        var_rgswbyrgsw_a += (
                RR(1 << (self.rgsw_by_rgsw_decomposer.ignore_bits_a()*2))/12
            *   var_sk_rlwe
            *   N
        )
        #print(f"RGSW x RGSW part A inexact noise std: {format_rr(sqrt(var_rgswbyrgsw_a-tmp))}")
        var_rgswbyrgsw_b = (d_b_rgsw_by_rgsw * ((B_rgsw_rgsw*B_rgsw_rgsw)/RR(12) * var_fresh * N))
        var_rgswbyrgsw_b += (
                RR(1 << (self.rgsw_by_rgsw_decomposer.ignore_bits_b()*2))/12
        )
        var_brk = k*(var_rgswbyrgsw_a+var_rgswbyrgsw_b)


        # RLWE x RGSW where RGSW has var_brk error variance
        d_a_rlwe_by_rgsw = self.rlwe_by_rgsw_decomposer.d_a
        d_b_rlwe_by_rgsw = self.rlwe_by_rgsw_decomposer.d_b
        B_rlwe_rgsw = RR(1<<self.rlwe_by_rgsw_decomposer.logB)
        var_rlwe_by_rgsw_a = (d_a_rlwe_by_rgsw * ((B_rlwe_rgsw*B_rlwe_rgsw)/12) * (var_brk) * N) 
        tmp = var_rlwe_by_rgsw_a
        #print(f"RLWE x RGSW Part A ks noise std: {format_rr(sqrt(tmp))}")
        var_rlwe_by_rgsw_a += (
                RR(1 << (self.rlwe_by_rgsw_decomposer.ignore_bits_a()*2))/12
            *   var_sk_rlwe
            *   N
        )
        #print(f"RLWE x RGSW Part A inexact noise std: {format_rr(sqrt(var_rlwe_by_rgsw_a-tmp))}")
        var_rlwe_by_rgsw_b = (d_b_rlwe_by_rgsw * ((B_rgsw_rgsw*B_rgsw_rgsw)/12) * (var_brk) * N) 
        var_rlwe_by_rgsw_b += (
                RR(1 << (self.rlwe_by_rgsw_decomposer.ignore_bits_b()*2))/12
        )
        var_rlwe_by_rgsw = (var_rlwe_by_rgsw_a+var_rlwe_by_rgsw_b)
        # print("var_rlwe_by_rgsw_a: ", sqrt(var_rlwe_by_rgsw_a))
        # print("var_rlwe_by_rgsw_b: ", sqrt(var_rlwe_by_rgsw_b))


        # Auto
        B_auto = RR(1<<self.auto_decomposer.logB)
        d_auto = self.auto_decomposer.d_a
        var_auto = ((B_auto*B_auto)/12)*(k*var)*N*d_auto
        var_auto += N*RR(1<<(self.auto_decomposer.ignore_bits_a()*2))/12*var_sk_rlwe


        # LWE ksk from rlwe secret to lwe secret
        B_lwe = RR(1<<self.lwe_decomposer.logB)
        d_lwe  = self.lwe_decomposer.d_a
        var_ks = (((B_lwe*B_lwe)/12)*(k*var)*d_lwe*N) 
        tmp = var_ks
        #print(f"LWE ks noise std: {format_rr(sqrt(tmp))}")
        var_ks += N*(RR(1<<(self.lwe_decomposer.ignore_bits_a()*2))/12*var_sk_rlwe)
        #print(f"LWE inexact noise std: {format_rr(sqrt(var_ks-tmp))}")


        # var ms1: Q -> Q_ks
        var_ms1 = ((N*var_sk_rlwe)+1) * RR(1/12)

        # var ms2: Q_ks -> q (odd mod switch. Set variance of rounding error 4/12)
        var_ms2 = ((n*var_sk_lwe)+1) * RR(4/12)

        Q_sq = RR(self.Q*self.Q)
        Q_ks_sq = RR(self.Q_ks*self.Q_ks)
        q_sq = RR(self.q*self.q)
        
        w = RR(self.w)

        msg = RR(pow(2, self.msg))
        padding = RR(pow(2, self.padding))
        
        # var_acc = (n*var_rlwe_by_rgsw)+(self.q*var_auto)
        worst_case_autos = (((w - 1)/w)*n)+(1/w)*(self.N)
        var_acc = n*var_rlwe_by_rgsw+var_auto*worst_case_autos
    
        #print(format_rr(sqrt((q_sq*(2*var_acc))/Q_sq)), format_rr(sqrt((q_sq*(var_ms1+var_ks))/Q_ks_sq)), format_rr(sqrt(var_ms1)), format_rr(sqrt(var_ms2)))

        var_zeta = (q_sq/Q_sq)*2*var_acc + (q_sq/Q_ks_sq)*(var_ms1+var_ks) + var_ms2 + RR(1/12)

        # [padding|msg|msg|0...0|noise]
        fail_prob = (RR(self.q)/(padding*msg**2)/sqrt(2*var_zeta)).erfc()
        """
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
            std_zeta: {format_rr(sqrt(var_zeta))}
            failure probability : {format_rr(fail_prob)}
        ''')
        """

        # if fail_prob != D(0):
        #print(f'message space: Z_{1<<self.msg}')
        #print(f'Failure probability log 2: {fail_prob.log2()}')

        return fail_prob.log2()

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