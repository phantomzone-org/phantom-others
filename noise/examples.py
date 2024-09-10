from parameters import Parameters, ParameterVariant
from decomposer import Decomposer
from secret import Secret
from math import log

logQ = 64
msg_bits = 1
lut_way = 3
padding = 1
N = 1<<13
max_log_B = 17
parties = 4
fresh_noise_std_lwe = 5.1
w = 20

P4_MSG1_LUT3_PADD1_80 = Parameters(
    msg_bits=msg_bits,
    lut_way=lut_way,
    padding=padding,
    logQ=logQ, 
    logQ_ks=25,
    logq=14, 
    lwe_sk=Secret.TernarySecret(N=590),
    rlwe_sk=Secret.TernarySecret(N),
    fresh_noise_std_lwe=fresh_noise_std_lwe,
    w=w,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=max_log_B,
    parties=parties,
)

print("P4_MSG1_LUT3_PADD1_80")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_80.optimize(-40)}')
P4_MSG1_LUT3_PADD1_80.print_decomposers()
#P4_MSG1_LUT3_PADD1_80.print_key_statistics()
print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_80.blind_rotate_ops_zp(), 2))
#P4_MSG1_LUT3_PADD1_80.security()
"""
arora-gb             :: rop: ≈2^inf, m: ≈2^430.4, dreg: 288, t: 22, mem: ≈2^inf, tag: arora-gb, ↻: ≈2^421.2, ζ: 297, |S|: ≈2^367.8, prop: ≈2^-51.2
bkw                  :: rop: ≈2^142.2, m: ≈2^129.1, mem: ≈2^130.1, b: 5, t1: 0, t2: 34, ℓ: 4, #cod: 536, #top: 0, #test: 54, tag: coded-bkw
usvp                 :: rop: ≈2^83.0, red: ≈2^83.0, δ: 1.006633, β: 183, d: 1116, tag: usvp
bdd                  :: rop: ≈2^80.6, red: ≈2^80.4, svp: ≈2^77.7, β: 174, η: 210, d: 1144, tag: bdd
dual                 :: rop: ≈2^84.8, mem: ≈2^46.8, m: 592, β: 189, d: 1182, ↻: 1, tag: dual
dual_hybrid          :: rop: ≈2^81.6, red: ≈2^81.5, guess: ≈2^75.3, β: 178, p: 4, ζ: 10, t: 20, β': 178, N: ≈2^35.3, m: 590
"""

P4_MSG1_LUT3_PADD1_100 = Parameters(
    msg_bits=msg_bits,
    lut_way=lut_way,
    padding=padding,
    logQ=logQ, 
    logQ_ks=25,
    logq=14, 
    lwe_sk=Secret.TernarySecret(N=735),
    rlwe_sk=Secret.TernarySecret(N),
    fresh_noise_std_lwe=fresh_noise_std_lwe,
    w=w,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=max_log_B,
    parties=parties,
)

print("P4_MSG1_LUT3_PADD1_100")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_100.optimize(-40)}')
P4_MSG1_LUT3_PADD1_100.print_decomposers()
#P4_MSG1_LUT3_PADD1_100.print_key_statistics()
print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_100.blind_rotate_ops_zp(), 2))
#P4_MSG1_LUT3_PADD1_100.security()
"""
arora-gb             :: rop: ≈2^inf, m: ≈2^533.3, dreg: 354, t: 22, mem: ≈2^inf, tag: arora-gb, ↻: ≈2^524.1, ζ: 370, |S|: ≈2^458.3, prop: ≈2^-63.6
bkw                  :: rop: ≈2^167.9, m: ≈2^154.3, mem: ≈2^155.3, b: 6, t1: 0, t2: 39, ℓ: 5, #cod: 670, #top: 0, #test: 65, tag: coded-bkw
usvp                 :: rop: ≈2^103.4, red: ≈2^103.4, δ: 1.005429, β: 252, d: 1383, tag: usvp
bdd                  :: rop: ≈2^100.9, red: ≈2^100.8, svp: ≈2^97.6, β: 243, η: 278, d: 1391, tag: bdd
dual                 :: rop: ≈2^105.5, mem: ≈2^61.8, m: 722, β: 259, d: 1457, ↻: 1, tag: dual
dual_hybrid          :: rop: ≈2^100.6, red: ≈2^100.5, guess: ≈2^96.3, β: 242, p: 4, ζ: 20, t: 20, β': 242, N: ≈2^49.6, m: 735
"""

P4_MSG1_LUT3_PADD1_128 = Parameters(
    msg_bits=msg_bits,
    lut_way=lut_way,
    padding=padding,
    logQ=logQ, 
    logQ_ks=25,
    logq=14, 
    lwe_sk=Secret.TernarySecret(N=925),
    rlwe_sk=Secret.TernarySecret(N),
    fresh_noise_std_lwe=fresh_noise_std_lwe,
    w=w,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=max_log_B,
    parties=parties,
)

print("P4_MSG1_LUT3_PADD1_128")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_128.optimize(-40)}')
P4_MSG1_LUT3_PADD1_128.print_decomposers()
#P4_MSG1_LUT3_PADD1_128.print_key_statistics()
print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_128.blind_rotate_ops_zp(), 2))
#P4_MSG1_LUT3_PADD1_128.security()

"""
arora-gb             :: rop: ≈2^inf, m: ≈2^668.5, dreg: 437, t: 23, mem: ≈2^inf, tag: arora-gb, ↻: ≈2^658.0, ζ: 464, |S|: ≈2^577.2, prop: ≈2^-78.6
bkw                  :: rop: ≈2^217.5, m: ≈2^203.9, mem: ≈2^204.9, b: 8, t1: 0, t2: 30, ℓ: 7, #cod: 829, #top: 0, #test: 96, tag: coded-bkw
usvp                 :: rop: ≈2^131.2, red: ≈2^131.2, δ: 1.004399, β: 346, d: 1737, tag: usvp
bdd                  :: rop: ≈2^128.5, red: ≈2^128.3, svp: ≈2^125.3, β: 336, η: 373, d: 1752, tag: bdd
dual                 :: rop: ≈2^133.9, mem: ≈2^82.1, m: 890, β: 355, d: 1815, ↻: 1, tag: dual
dual_hybrid          :: rop: ≈2^126.9, red: ≈2^126.9, guess: ≈2^117.7, β: 331, p: 3, ζ: 20, t: 40, β': 331, N: ≈2^67.8, m: 925
"""