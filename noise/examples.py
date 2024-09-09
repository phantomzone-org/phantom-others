from parameters import Parameters, ParameterVariant
from decomposer import Decomposer
from secret import Secret
from math import log

P2_MSG2_PADD1 = Parameters(
    msg_bits=2,
    padding=1,
    logQ=64, 
    logQ_ks=32,
    logq=14, 
    lwe_sk=Secret.ErrorDistribution(N=600),
    rlwe_sk=Secret.TernarySecret(N=1<<12),
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=17,
    parties=4,
)

print(f'PR[e > Q/2^{{2*msg+padding}} | drift > N/2^{{msg}}]=2^{P2_MSG2_PADD1.optimize(-80)}')
P2_MSG2_PADD1.print_decomposers()
P2_MSG2_PADD1.print_key_statistics()
print("log2(#ops_Zp):", log(P2_MSG2_PADD1.blind_rotate_ops_zp(), 2))
#P2_MSG2_PADD1.security()