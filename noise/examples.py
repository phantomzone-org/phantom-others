from parameters import Parameters, ParameterVariant
from decomposer import Decomposer
from secret import Secret

P2_MSG2_PADD1 = Parameters(
    msg=2,
    padding=1,
    logQ=54, 
    logQ_ks=18,
    logq=16, 
    lwe_sk=Secret.ErrorDistribution(N=590),
    rlwe_sk=Secret.TernarySecret(N=1<<11),
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    parties=2,
)

print(f'PR[fail]=2^{P2_MSG2_PADD1.optimize(-128)}')
P2_MSG2_PADD1.print_decomposers()
