from parameters import Parameters, ParameterVariant
from decomposer import Decomposer
from secret import Secret
from math import log

P2_MSG2_PADD1 = Parameters(
    msg=2,
    padding=1,
    logQ=54, 
    logQ_ks=18,
    logq=13, 
    lwe_sk=Secret.ErrorDistribution(N=590),
    rlwe_sk=Secret.TernarySecret(N=1<<12),
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=2,
)

print(f'PR[e > Q/2^{{2*msg+padding}} | drift > N/2^{{msg}}]=2^{P2_MSG2_PADD1.optimize(-128)}')
P2_MSG2_PADD1.print_decomposers()
print("log2(#ops_Zp):", log(P2_MSG2_PADD1.blind_rotate_ops_zp(), 2))
#P2_MSG2_PADD1.security()