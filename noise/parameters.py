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
        msg_bits,
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
        max_logB = None,
        rgsw_by_rgsw_decomposer: Decomposer=None,
        rlwe_by_rgsw_decomposer: Decomposer=None,
        auto_decomposer: Decomposer=None,
        lwe_decomposer:Decomposer=None,
        non_interactive_uitos_decomposer:Decomposer=None,
    ):
        self.k = parties
        self.msg_bits = msg_bits
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
                max_logB=max_logB,
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
                max_logB=max_logB,
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

        self.Q = RR(1<<logQ)
        self.Q_ks = RR(1<<logQ_ks)
        self.q = RR(1<<logq)

        self.N = RR(rlwe_sk.dimension)
        self.n = RR(lwe_sk.dimension)

        self.fresh_noise_std = RR(fresh_noise_std)
        self.fresh_noise_var = self.fresh_noise_std*self.fresh_noise_std

        self.w = RR(w)

        self.variant = variant

        self.lwe_sk = lwe_sk
        self.rlwe_sk = rlwe_sk

        self.var_sk_rlwe = RR(self.k)*self.rlwe_sk.variance()
        self.var_sk_lwe = RR(self.k)*self.lwe_sk.variance()

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

        # In order of affecting the complexity of the bootstrapping
        pr_fail = self.rlwe_by_rgsw_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.auto_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.lwe_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.rgsw_by_rgsw_decomposer.optimize(self.pr_fail, log_fail)
        
        # Second pass
        pr_fail = self.rlwe_by_rgsw_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.auto_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.lwe_decomposer.optimize(self.pr_fail, log_fail)
        pr_fail = self.rgsw_by_rgsw_decomposer.optimize(self.pr_fail, log_fail)
        
        return pr_fail

    def print_decomposers(self):
        print(">>> DECOMPOSERS <<<")
        print(f'rgsw_by_rgsw_decomposer: logB={self.rgsw_by_rgsw_decomposer.logB:<2d} d_a={self.rgsw_by_rgsw_decomposer.d_a:<2d} d_b={self.rgsw_by_rgsw_decomposer.d_b:<2d}')
        print(f'rlwe_by_rgsw_decomposer: logB={self.rlwe_by_rgsw_decomposer.logB:<2d} d_a={self.rlwe_by_rgsw_decomposer.d_a:<2d} d_b={self.rlwe_by_rgsw_decomposer.d_b:<2d}')
        print(f'auto_decomposer        : logB={self.auto_decomposer.logB:<2d} d_a={self.auto_decomposer.d_a:<2d}')
        print(f'lwe_decomposer         : logB={self.lwe_decomposer.logB:<2d} d_a={self.lwe_decomposer.d_a:<2d}')

    def print_key_statistics(self):
        print(">>> KEY STATISTICS <<<")
        print("RGSW log2 var:", log(self.var_brk(), 2))
        print("rlwe x auto log2 var:", log(self.var_auto(), 2))
        print("rlwe x RGSW log2 var", log(self.var_rlwe_rgsw(self.var_brk()), 2))

    def var_rgsw_fresh(self):
        n = RR(self.n)
        N = RR(self.N)
        k = RR(self.k)

        var_fresh = RR(0)
        match self.variant:
            case ParameterVariant.INTERACTIVE_MULTIPARTY:
                var_fresh = self.fresh_noise_var * (2 * self.N * self.var_sk_rlwe+1)
            case ParameterVariant.NON_INTERACTIVE_MULTIPARTY:
                B = RR(1<<self.non_interactive_uitos_decomposer.logB)
                var_fresh = self.N * self.k * self.fresh_noise_var * B**2 / RR(12) * RR(self.non_interactive_uitos_decomposer.d_a)
                var_fresh += self.N * self.var_sk_rlwe * RR(1<<(self.non_interactive_uitos_decomposer.ignore_bits_a()*2))/RR(12)
                var_fresh += self.fresh_noise_var * self.var_sk_rlwe
        return var_fresh

    def var_rgsw_rgsw(self):

        var_rgsw_fresh = self.var_rgsw_fresh()

        B_rgsw_rgsw = RR(1<<self.rgsw_by_rgsw_decomposer.logB)

        # Approximation error induced by ignoring some least signifcant bits. 
        # The variance of ignored bits is (2^{ignored_bits})^2
        var_rgswbyrgsw_a  = self.N * var_rgsw_fresh * self.rgsw_by_rgsw_decomposer.d_a * (B_rgsw_rgsw**2)/RR(12)
        var_rgswbyrgsw_b  = self.N * var_rgsw_fresh * self.rgsw_by_rgsw_decomposer.d_b * (B_rgsw_rgsw**2)/RR(12)

        var_rgswbyrgsw_a += self.N * self.var_sk_rlwe * RR(1 << (self.rgsw_by_rgsw_decomposer.ignore_bits_a()*2))/RR(12)
        var_rgswbyrgsw_b += RR(1 << (self.rgsw_by_rgsw_decomposer.ignore_bits_b()*2))/RR(12)

        return var_rgswbyrgsw_a + var_rgswbyrgsw_b

    def var_brk(self):
        return self.var_rgsw_fresh() + (self.k-1) * self.var_rgsw_rgsw()

    def var_rlwe_rgsw(self, var_rgsw):
        B_rlwe_rgsw = RR(1<<self.rlwe_by_rgsw_decomposer.logB)

        var_brk = self.var_brk()

        var_rlwe_by_rgsw_a = self.N * var_brk * self.rlwe_by_rgsw_decomposer.d_a * (B_rlwe_rgsw**2)/RR(12)
        var_rlwe_by_rgsw_b = self.N * var_brk * self.rlwe_by_rgsw_decomposer.d_b * (B_rlwe_rgsw**2)/RR(12) 

        var_rlwe_by_rgsw_a += self.N * self.var_sk_rlwe * RR(1 << (self.rlwe_by_rgsw_decomposer.ignore_bits_a()*2))/RR(12)
        var_rlwe_by_rgsw_b += RR(1 << (self.rlwe_by_rgsw_decomposer.ignore_bits_b()*2))/RR(12)

        return var_rlwe_by_rgsw_a + var_rlwe_by_rgsw_b

    def var_acc(self):
        return self.n * self.var_rlwe_rgsw(self.var_brk()) + self.worst_case_autos() * self.var_auto()

    def var_auto(self):
        var_auto  = self.N * self.k * self.auto_decomposer.d_a * RR(1<<(self.auto_decomposer.logB*2))/RR(12)
        var_auto += self.N * self.var_sk_rlwe * RR(1<<(self.auto_decomposer.ignore_bits_a()*2))/RR(12)
        return  var_auto

    def worst_case_autos(self):
        return (((self.w - 1)/self.w)*self.n)+(1/self.w)*(self.N)

    def var_ks_rlwe_to_lwe(self):
        var_ks  = self.N * self.k * self.fresh_noise_var *self.lwe_decomposer.d_a * RR(1<<(self.auto_decomposer.logB*2))/RR(12)
        var_ks += self.N * self.var_sk_rlwe * RR(1<<(self.lwe_decomposer.ignore_bits_a()*2))/RR(12)
        return var_ks

    def var_blind_rotate(self):

        Q_sq = self.Q*self.Q
        Q_ks_sq = self.Q_ks*self.Q_ks
        q_sq = self.q*self.q

        var_ms1 = ((self.N*self.var_sk_rlwe)+1) * RR(1/12)  # Q -> Q_ks
        var_ms2 = ((self.n*self.var_sk_lwe)+1) * RR(1/12) # Q_ks -> q
        
        # blind_rotate(ks_rlwe_to_lwe(ct0 + msg * ct1)) + [Q -> Q_ks] + [Q_ks -> q]
        return (q_sq/Q_sq)*self.var_acc()*(1+pow(2, self.msg_bits)) + (q_sq/Q_ks_sq)*(var_ms1+self.var_ks_rlwe_to_lwe()) + var_ms2 + RR(1/12)

    def auto_ops_zp(self):
        return self.N * self.auto_decomposer.d_a * (log(self.N, 2) + 1)

    def rlwe_rgsw_ops_zp(self):
        return self.N * (self.rlwe_by_rgsw_decomposer.d_a+self.rlwe_by_rgsw_decomposer.d_b) * (log(self.N, 2)+2)

    def blind_rotate_ops_zp(self):
        return self.n * self.rlwe_rgsw_ops_zp() + self.worst_case_autos() * self.auto_ops_zp()

    def var_drift(self):
        return RR(1/12) + (self.n*self.var_sk_lwe) * RR(4/12) #  odd mod switch

    def blind_rotate_fail_prob(self):
        msg_space = pow(2, self.padding)*pow(4, self.msg_bits) # [padding|msg|msg]
        return ((RR(self.q)/msg_space)/sqrt(2*self.var_blind_rotate())).erfc() # PR[e > q/msg_space]

    def drift_fail_prob(self):
        return ((self.N/pow(2, self.msg_bits))/sqrt(2*self.var_drift())).erfc() # PR[e > N/msg^2]

    def pr_fail(self):
        return (1-(1-self.blind_rotate_fail_prob())*(1-self.drift_fail_prob())).log2()

    def security(self):
        # LWE
        lwe = LWE.Parameters(n=int(self.n), q=int(self.Q_ks), Xs=self.lwe_sk.distr, Xe=ND.DiscreteGaussian(self.fresh_noise_std), m=int(2*self.n))
        lwe_res = LWE.estimate(lwe, red_cost_model = RC.BDGL16)

        print("LWE Security")
        print(lwe)
        print(lwe_res)
        
        print("")

        rlwe = LWE.Parameters(n=int(self.N), q=int(self.Q), Xs=self.rlwe_sk.distr, Xe=ND.DiscreteGaussian(self.fresh_noise_std), m=int(2*self.N))
        rlwe_res = LWE.estimate(rlwe, red_cost_model = RC.BDGL16)

        print("RLWE Security")
        print(rlwe)
        print(rlwe_res)
