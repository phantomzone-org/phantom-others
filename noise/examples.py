from parameters import Parameters, ParameterVariant
from decomposer import Decomposer
from secret import Secret
from math import log

# logQ = 64
# msg_bits = 1
# lut_way = 3
# padding = 1
# N = 1<<13
# max_log_B = 17
# parties = 4
# fresh_noise_std_lwe = 5.1
# w = 10

P4_BOOL_60 = Parameters(
    msg_bits=0,
    lut_way=0,
    padding=0,
    logQ=64, 
    logQ_ks=18,
    logq=10, 
    lwe_sk=Secret.TernarySecret(N=300),
    rlwe_sk=Secret.TernarySecret(1<<10),
    w=10,
    variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
    max_logB=18,
    # rgsw_by_rgsw_decomposer=Decomposer.double_decomposer(
    #     logQ=64, 
    #     logB=14,
    #     d_a=3, 
    #     d_b=3,
    # ),
    parties=4,
)

print("P4_BOOL_60")
print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_BOOL_60.optimize(-50)}')
P4_BOOL_60.print_decomposers()
#P4_MSG1_LUT3_PADD1_80.print_key_statistics()
print("log2(#ops_Zp):", log(P4_BOOL_60.blind_rotate_ops_zp(), 2))
print(f"Decryption failure probability: 2^{P4_BOOL_60.decryption_failure(pack_lwe=True, packing_logQ=61, logB=20, d=1)}")
P4_BOOL_60.security()

# P4_MSG1_LUT3_PADD1_80 = Parameters(
#     msg_bits=msg_bits,
#     lut_way=lut_way,
#     padding=padding,
#     logQ=logQ, 
#     logQ_ks=18,
#     logq=14, 
#     lwe_sk=Secret.TernarySecret(N=411),
#     rlwe_sk=Secret.TernarySecret(N),
#     fresh_noise_std_lwe=fresh_noise_std_lwe,
#     w=w,
#     variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
#     max_logB=max_log_B,
#     parties=parties,
# )

# print("P4_MSG1_LUT3_PADD1_80")
# print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_80.optimize(-40)}')
# P4_MSG1_LUT3_PADD1_80.print_decomposers()
# #P4_MSG1_LUT3_PADD1_80.print_key_statistics()
# print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_80.blind_rotate_ops_zp(), 2))
# #P4_MSG1_LUT3_PADD1_80.security()
# """
# arora-gb             :: rop: ≈2^inf, m: ≈2^304.1, dreg: 207, t: 22, mem: ≈2^810.6, tag: arora-gb, ↻: ≈2^294.9, ζ: 208, |S|: ≈2^256.0, prop: ≈2^-36.7
# bkw                  :: rop: ≈2^122.8, m: ≈2^111.1, mem: ≈2^112.1, b: 6, t1: 0, t2: 17, ℓ: 5, #cod: 358, #top: 0, #test: 54, tag: coded-bkw
# usvp                 :: rop: ≈2^83.1, red: ≈2^83.1, δ: 1.006589, β: 185, d: 800, tag: usvp
# bdd                  :: rop: ≈2^80.1, red: ≈2^79.8, svp: ≈2^78.0, β: 174, η: 211, d: 748, tag: bdd
# dual                 :: rop: ≈2^86.0, mem: ≈2^48.1, m: 413, β: 195, d: 824, ↻: 1, tag: dual
# dual_hybrid          :: rop: ≈2^81.0, red: ≈2^81.0, guess: ≈2^75.2, β: 178, p: 4, ζ: 10, t: 20, β': 178, N: ≈2^36.4, m: 411
# """

# P4_MSG1_LUT3_PADD1_100 = Parameters(
#     msg_bits=msg_bits,
#     lut_way=lut_way,
#     padding=padding,
#     logQ=logQ, 
#     logQ_ks=18,
#     logq=14, 
#     lwe_sk=Secret.TernarySecret(N=515),
#     rlwe_sk=Secret.TernarySecret(N),
#     fresh_noise_std_lwe=fresh_noise_std_lwe,
#     w=w,
#     variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
#     max_logB=max_log_B,
#     parties=parties,
# )

# print("P4_MSG1_LUT3_PADD1_100")
# print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_100.optimize(-40)}')
# P4_MSG1_LUT3_PADD1_100.print_decomposers()
# #P4_MSG1_LUT3_PADD1_100.print_key_statistics()
# print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_100.blind_rotate_ops_zp(), 2))
# #P4_MSG1_LUT3_PADD1_100.security()
# """
# arora-gb             :: rop: ≈2^inf, m: ≈2^377.7, dreg: 254, t: 22, mem: ≈2^1008.4, tag: arora-gb, ↻: ≈2^368.4, ζ: 260, |S|: ≈2^320.2, prop: ≈2^-46.0
# bkw                  :: rop: ≈2^141.7, m: ≈2^129.4, mem: ≈2^130.3, b: 7, t1: 0, t2: 20, ℓ: 6, #cod: 451, #top: 0, #test: 65, tag: coded-bkw
# usvp                 :: rop: ≈2^103.5, red: ≈2^103.5, δ: 1.005401, β: 254, d: 971, tag: usvp
# bdd                  :: rop: ≈2^100.2, red: ≈2^99.9, svp: ≈2^97.9, β: 242, η: 279, d: 936, tag: bdd
# dual                 :: rop: ≈2^106.8, mem: ≈2^63.0, m: 501, β: 265, d: 1016, ↻: 1, tag: dual
# dual_hybrid          :: rop: ≈2^99.8, red: ≈2^99.7, guess: ≈2^96.1, β: 241, p: 4, ζ: 20, t: 20, β': 241, N: ≈2^49.4, m: 515
# """

# P4_MSG1_LUT3_PADD1_128 = Parameters(
#     msg_bits=msg_bits,
#     lut_way=lut_way,
#     padding=padding,
#     logQ=logQ, 
#     logQ_ks=18,
#     logq=14, 
#     lwe_sk=Secret.TernarySecret(N=655),
#     rlwe_sk=Secret.TernarySecret(N),
#     fresh_noise_std_lwe=fresh_noise_std_lwe,
#     w=w,
#     variant=ParameterVariant.INTERACTIVE_MULTIPARTY,
#     max_logB=max_log_B,
#     parties=parties,
# )

# print("P4_MSG1_LUT3_PADD1_128")
# print(f'PR[e > Q/{{msg_space}} | drift > N/{{2*msg_space}}]=2^{P4_MSG1_LUT3_PADD1_128.optimize(-40)}')
# P4_MSG1_LUT3_PADD1_128.print_decomposers()
# #P4_MSG1_LUT3_PADD1_128.print_key_statistics()
# print("log2(#ops_Zp):", log(P4_MSG1_LUT3_PADD1_128.blind_rotate_ops_zp(), 2))
# P4_MSG1_LUT3_PADD1_128.security()

# """
# arora-gb             :: rop: ≈2^inf, m: ≈2^476.7, dreg: 317, t: 22, mem: ≈2^inf, tag: arora-gb, ↻: ≈2^467.5, ζ: 330, |S|: ≈2^407.0, prop: ≈2^-58.3
# bkw                  :: rop: ≈2^177.9, m: ≈2^165.3, mem: ≈2^166.3, b: 9, t1: 0, t2: 20, ℓ: 8, #cod: 576, #top: 0, #test: 84, tag: coded-bkw
# usvp                 :: rop: ≈2^131.8, red: ≈2^131.8, δ: 1.004365, β: 350, d: 1208, tag: usvp
# bdd                  :: rop: ≈2^128.5, red: ≈2^128.3, svp: ≈2^125.3, β: 338, η: 373, d: 1179, tag: bdd
# dual                 :: rop: ≈2^135.4, mem: ≈2^83.6, m: 617, β: 362, d: 1272, ↻: 1, tag: dual
# dual_hybrid          :: rop: ≈2^126.1, red: ≈2^126.0, guess: ≈2^117.6, β: 330, p: 3, ζ: 20, t: 40, β': 330, N: ≈2^67.7, m: 655
# """