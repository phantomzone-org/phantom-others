from parameters import Parameters, ParameterVariant
from decomposer import Decomposer
from secret import Secret
from math import log


P4_MSG1_LUT3_PADD1_80 = Parameters(
    msg_bits=1,
    lut_way=3,
    padding=1,
    logQ=54, 
    logQ_ks=20,
    logq=14, 
    lwe_sk=Secret.TernarySecret(N=480),
    rlwe_sk=Secret.TernarySecret(N=1<<13),
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=17,
    parties=4,
)

print("P4_MSG1_LUT3_PADD1_80")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_80.optimize(-40)}')
P4_MSG1_LUT3_PADD1_80.print_decomposers()
#P4_MSG1_LUT3_PADD1_80.print_key_statistics()
print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_80.blind_rotate_ops_zp(), 2))
#P4_MSG1_LUT3_PADD1_80.security()

P4_MSG1_LUT3_PADD1_100 = Parameters(
    msg_bits=1,
    lut_way=3,
    padding=1,
    logQ=54, 
    logQ_ks=20,
    logq=14, 
    lwe_sk=Secret.TernarySecret(N=608),
    rlwe_sk=Secret.TernarySecret(N=1<<13),
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=17,
    parties=4,
)

print("P4_MSG1_LUT3_PADD1_100")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_100.optimize(-40)}')
P4_MSG1_LUT3_PADD1_100.print_decomposers()
#P4_MSG1_LUT3_PADD1_100.print_key_statistics()
print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_100.blind_rotate_ops_zp(), 2))
#P4_MSG1_LUT3_PADD1_100.security()

P4_MSG1_LUT3_PADD1_128 = Parameters(
    msg_bits=1,
    lut_way=3,
    padding=1,
    logQ=54, 
    logQ_ks=24,
    logq=14, 
    lwe_sk=Secret.TernarySecret(N=912),
    rlwe_sk=Secret.TernarySecret(N=1<<13),
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=17,
    parties=4,
)

print("P4_MSG1_LUT3_PADD1_128")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_128.optimize(-40)}')
P4_MSG1_LUT3_PADD1_128.print_decomposers()
#P4_MSG1_LUT3_PADD1_128.print_key_statistics()
print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_128.blind_rotate_ops_zp(), 2))
#P4_MSG1_LUT3_PADD1_128.security()

