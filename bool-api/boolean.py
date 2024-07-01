from __future__ import annotations
import numpy as np
import copy
import math

def arbitrary_bit_equality(a: [bool], b: [bool]) -> bool:
    '''
    Returns True if a == b, False otherwise. 

    Requires: 
        - N XNORs
        - N-1 ANDs
    where N is bit width
    '''
    assert len(a) == len(b)
    out = not (a[0]^b[0])
    for i in range(1, len(a)):
        # XNOR a[i], b[i]
        out = out & (not (a[i]^b[i]))
    return out

def arbitrary_signed_bit_comparator(a: [bool], b:[bool]) -> bool:
    a = copy.deepcopy(a)
    b = copy.deepcopy(b)
    a[-1]  = not a[-1]
    b[-1]  = not b[-1]

    return arbitrary_unsigned_bit_comparator(a=a, b=b)

def arbitrary_unsigned_bit_comparator(a: [bool], b:[bool]) -> bool:
    '''
    Assumes A and B are unsigned and returns True if A > B, otherwise False

    Uses two gates: 
    (1) a[i] & !(b[i])
        Outputs 1 iif a[i] > b[i]. In other words, outputs 1 only when a[i] = 1 and b[i] = 0

    (2) !(a[i]^b[i])
        Outputs 1 iif a[i] == b[i], otherwise 0

    We then combine gates (1) and (2) to output whether A > B. The trick is to evaluate (1) for each 
    bit index, but to consider output for bit i valid only when bits at higher signicant positions are 
    equal. This follows from the simple fact that first non-differing bit is enough to decide whether A>B. 
    Thus we simply ignore output of (1) bits at lesser signifcanct position the moment we find a position 
    with unequal bits

    This function can be used to evaluate A<B, A<=B, A>=B
    - A<B: Call the function with reversed position
    - A<=B: !(A>B)
    - A>=B: !(A<B)
    
    # Signed comparison

    Signed comparison is equivalent to Unsigned comparison after flipping the MSB. In other words, order of signed 
    values with flipped MSB equals order of unsigned values. 

    For ex,

    GIVE 8 BIT EXAMPLE HERE
    '''
    N = len(a)
    assert len(a) == len(b)

    # N-1
    comp_bit = a[N-1] & (not b[N-1]) # N-1

    # N-2
    casc_bit = not (a[N-1]^b[N-1])
    comp_bit = comp_bit | ((a[N-2] & (not b[N-2])) & casc_bit)

    for j in range(N-3, -1, -1):
        casc_bit = casc_bit & (not (a[j+1]^b[j+1]))
        comp_bit = comp_bit | ((a[j] & (not b[j])) & casc_bit)

    return comp_bit

def half_adder(A: bool, B: bool) -> (bool, bool):
    '''
    Adds two bits A and B and returns sum S and Carry C. 
    '''

    out = A ^ B
    carry = A & B
    
    return (out, carry)

def full_adder(A: bool, B: bool, carry_in: bool) -> (bool, bool):
    '''
    Full adder. Returns (Sum, Carry out)
    '''

    # A xor B
    A_xor_B = A ^ B

    # S = (A xor B) xor C
    out = A_xor_B ^ carry_in

    # A and B
    A_and_B = A & B

    # (A xor B) and C
    A_xor_B__and_C = A_xor_B & carry_in

    # C_out = (A and B) or ((A xor B) and C)
    carry_out = A_and_B | A_xor_B__and_C
    
    return (out, carry_out)

def arbitrary_bit_adder(a: [bool], b: [bool], carry_in: bool) -> ([bool], bool, bool):
    '''
    Returns (sum = (a + b + c_in) mod 2^N, c_{N-1}, c_{N-2}) where N indicates no. of bits. 

    - c_{N-1}, c_{N-2} are returned to set necessary overflow flag in higher level APIs
    '''

    N = len(a)
    assert len(a) == len(b), f'len(a)={len(a)} != len(b)={len(b)}'
    assert N > 2
    
    out = [False for i in range(N)]

    start = 0
    carry = carry_in
    if carry_in == False: 
        # LSB has no carry_in
        (out[0], carry) = half_adder(a[0], b[0])
        start=1

    for i in range(start, N-2):
        (out[i], carry) = full_adder(a[i], b[i], carry_in=carry)   

    # handle 6,7 
    (out[N-2], c_Nminus2) = full_adder(a[N-2], b[N-2], carry_in=carry)   
    (out[N-1], c_Nminus1) = full_adder(a[N-1], b[N-1], carry_in=c_Nminus2)   

    return (out, c_Nminus1, c_Nminus2)

def arbitrary_bit_subtractor(a: [bool], b: [bool], borrow_in: bool) -> ([bool], bool, bool):
    '''
    Returns (a - b - borrow_in, c_{N-1}, c_{N-2}).

    - c_{N-1}, c_{N-2} are returned to set overflow flags in higher level APIs

    To calculate `a - b - borrow_in` observe: 
        a - b - borrow_in (mod 2^N)
        a + 2^N - b - borrow_in (mod 2^N)                         (2N - b: 2's complement of b)
        a + (2^N-1 - b) + 1 - borrow_in (mod 2^N)                 (2N - 1 - b: 1's complement of b)

    Hence, we calculate:
        a + (2^N-1 - b) + 1 - borrow_in (mod 2^N)
    and use arbitrary_bit_adder as a subroutine with carry_in flag set to 1 - borrow_in (i.e 1^borrow_in)

    Note that (2^N-1 - b) = !b 
    '''
    invert_b = [not i for i in b]
    carry_in = True^borrow_in
    return arbitrary_bit_adder(a=a, b=invert_b, carry_in=carry_in)

def eight_bits_mul(a: [bool], b: [bool]) -> [bool]:
    '''
    Returns a x b (mod 2^N).

    We implement vanilla decimal type long multiplication algorithm but for binary values. However, since the resulting value will be 2*(2^N) we avoid
    calculation of 2^N MSBs. 

    The algorithm when represented as a circuit is simple pattern with repeated invocations of half adders and full adders like a waterfall. Carries from last bits are
    used by the next. The function that I have written below implements that exactly. 

    Note that we delibrerately avoid calculating 2N MSBs because then no. of gates requires will be twice of what we have now. This implies we also lose the abiliyt to 
    determine, irrpesctive of signed/unsigned values, whether the output has overflown. When values are unsigned output overflows when result is > 2^N. When values are
    signed output overflows when result is not in range [-2^{N-1}, 2^{N-1}-1]. Overflow detection is left for future. 

    # References
    - https://inst.eecs.berkeley.edu/~eecs151/sp18/files/Lecture21.pdf
    - https://pages.cs.wisc.edu/%7Emarkhill/cs354/Fall2008/beyond354/int.mult.html
    '''
    assert len(a) == len(b) == 8

    sum = [False for i in range(8)]
    carries = [False for i in range(8)]
    for i in range(8):
        if i == 0:
            sum[0] = a[0]&b[0]
        elif i == 1:
            sum[1], carries[1] = half_adder(A=(a[1]&b[0]), B=(a[0]&b[1]))
        else:
            (tmp, carries[1]) = full_adder(A=(a[i]&b[0]), B=(a[i-1]&b[1]),carry_in=carries[1])
            for j in range(2, i):
                (tmp, carries[j]) = full_adder(A=tmp, B=(a[i-j]&b[j]),carry_in=carries[j])
            # half adder
            (sum[i], carries[i]) = half_adder(A=tmp, B=(a[0]&b[i]))

    return sum 

def absolute(a:[bool]) -> [bool]:
    '''
    Assumes `a` is in signed representation as 2s complement. Returns abs(a)

    Note that +ve counterpart of -ve value in 2s complement equals 2s complement. For example, 
    assume a is -ve. It is represented in 2s complement as:
        2^N - a (mod 2^N)
    Its +ve counterpart equals:
        2^N - (2^N - a) (mod 2^N) = a (mod 2^N)
    That is, 2s complement of (2^N - a)
    '''
    assert len(a) == 8

    # if a is negative then send it to its 2's complement (ie its +ve counterpart)
    a_if_neg = [not i for i in a] # 1's complement
    carry = True # +1 to send 1s complement to 2s complement
    for i in range(8):
        (a_if_neg[i], carry) = half_adder(A=a_if_neg[i], B=carry)

    # if a is -ve (the case when MSB of a is set to True), then return a_if_neg. Otherwise returns `a` as it is
    return mux_bool_vec(bit=a[-1], a=a_if_neg, b=a)


def mux_bool_vec(bit: bool, a:[bool], b:[bool]) -> [bool]:
    '''
    Muxer that returns `a` when bit=True, otherwise returns `b`
    '''
    if bit: 
        return copy.deepcopy(a)
    else:
        return copy.deepcopy(b)

def mux_bool(bit: bool, a:bool, b:bool) -> bool:
    '''
    Muxer that returns `a` when bit=True, otherwise returns `b`
    '''
    if bit: 
        return a
    else:
        return b

def is_zero(a: [bool]) -> bool:
    '''
    Returns True if a == 0, otherwise returns False

    Requires: 
        - N/2 + (N/2)-1 ANDs
    '''
    N = len(a)

    assert N & (N - 1)  == 0 

    out = ((not a[0]) & (not (a[(N>>1)])))
    for i in range(1, N>>1):
        out = out & ((not a[i]) & (not (a[i+(N>>1)])))

    return out 

def arbitrary_unsigned_division(a: [bool], b:[bool]) -> ([bool], [bool]):
    '''
    Returns c, d s.t.
        a = b*c + d
    That is, given dividend `a` divisor `b` returns quotient `c` and remainder `d`.

    Below we implement restoring division (based on long division of binary digits) for unsigned integers. 
    For signed integers, pass them as their absolute values to the function and adjust the signs of quotient and remainder
    as per dividiend and divisor. 

    There's a variant called non-restoring division that observes: given after each iteration, if remainder < divisor
    we restore the remainder and left shift. Otherwise we left shift `remainder-divisor`. That is in each iteration we do, 
        PrevR = 2R
        R = PrevR - D 
        if R < 0:
            R = PrevR
        else:
            R = R

    We can change this to:
        Initialise R = 2R - D then for each iteration, 
            if R < 0:
                R = 2R + D (2*(R+D)-D = 2R+D)
            else:
                R = 2R - D
    Apparently latter is more hardware friendly because there's no re-storing step. Hence, no need to reserve registers to restore PrevR.
    
    But the same does not apply for us. Non-restoring will require 1 addition and 1 subtraction per iteration (since we execture all branches and then mux
    based on whether R < 0). More so, detecting R < 0 in restoring division is as simple as checking overflow in unsigned subtraction. I suspect in non-restoring division,
    to check whether R < 0 one will have extend the bitwidth of R and D to 2N to prevent overflow for corner cases for R and D (since they are doubled every iteration). 
    Due to this in non-restoring divison every arithmetic operation is for bit-width 2N. Whereas in restoring division we can simply stick to N bit-width arithemtic
    operations. 

    # Div by zero
    When b=0, we call it div by zero. The calling function must indenepdently check whether b=0. When b=0, quotient=MAX value N bits can support and remainder=a.

    - Both `a` and `b` are assumed to be +ve
    - Division by 0 (i.e. b = 0) must be indepenedetly checked by the caller. 
    '''
    N = len(a)
    assert len(b) == N

    quotient = [False for i in range(8)]
    remainder = [False for i in range(8)]

    for i in range(N):
        # Like long division, we iterate from MSB to LSB

        remainder = [a[N-1-i]] + remainder[:N-1]
        # Unsigned subtraction remainder - divisor
        (check, c_last, c_slast) = arbitrary_bit_subtractor(a=remainder, b=b, borrow_in=False)
        # Unsigned subtraction overflows when c_{N-1} is not set. Overflows means remainder < divisor
        
        # We implement restoring long division. We require to roll back, that is leave remainder untouched, when remainder < divisor. 
        # We do this using a muxer. We select `remainder - divisor` (i.e. `check`) if the subtraction has not overflown (i.e c_last = True), otherwise we select
        # remainder (i.e. c_last=False)
        remainder = mux_bool_vec(bit=c_last, a=check, b=remainder)

        # After each iteration we need to update quotient bit. q[N-1-i] must be set to zero if roll back happnes otherwise q[N-1-i] = 1
        quotient[N-1-i] = mux_bool(bit=c_last, a=True, b=False)

    return (quotient, remainder)


class FheInt8:
    def __init__(self, bits: [bool]):
        assert len(bits) == 8
        self.bits = bits

    def from_int8(v: int) -> FheInt8:
        # value can be in range [-128, 127]
        if v > 0:
            assert v <= 127
        else:
            assert v >= -128       
        
        if v < 0:
            # send to 2's complement
            # invert bits
            v = ~abs(v)
            v += 1

        # extract 8 bits
        bits = []
        for i in range(0, 8):
            bits.append(((v >> i) & 1)==1)

        return FheInt8(bits=bits)
        
    def to_int8(self) -> int:
        bits = self.bits
        assert len(bits) == 8

        value = 0
        for i in range(0, 8):
            if bits[i] == True:
                value += (1 << i)
        
        if value > 127:
            return -(256 - value)
        else:
            return value

    def Add(self, b: FheInt8) -> (FheInt8, bool): 
        '''
        Adds two signed integers

        # Overflow
        Overflows only when the two signed integers added have same sign and sign of the ouput does not equals sign of the inputs. 
        This happens when c_7 XOR c_6 = 1, that is when either of them, but not both, are set to 1. 

        When adding two positive values, if c_6 is 1 then it means the addition has overflown into the -ve region. And if c_7 is 0, it means 
        the overflown value has not wrapped around back into +ve region. Hence, an overflow. 

        When adding two negative values, if c_6 is 0, and given that MSB of both addends is 1, MSB of output will be 0, meaning output has
        overflown, and c_7 = 1. If c_6 = 1, MSB of output will be 1 with c_7 = 1. 
        '''
        (out_bits, c_7, c_6) = arbitrary_bit_adder(a=self.bits, b=b.bits, carry_in=False)
        overflow = c_7^c_6
        return (FheInt8(bits=out_bits), overflow)

    def Sub(self, b: FheInt8) -> (FheInt8, bool):
        '''
        Subtracts two signed integers

        # Overflow
        Subtraction only overflows when the two inputs have opposite signs. Overflow conditions are same as for addition if one negates + values
        as 2s complement. 

        When subtracting -y from +x it equals addition of +x and +y
            x - (-y) => x + y

        When subtracting +y fron -x it equals addition of -x and -y
            -x - (+y) => - x - y

        '''
        (out_bits, c_7, c_6) = arbitrary_bit_subtractor(a=self.bits, b=b.bits, borrow_in=False)
        overflow = c_7^c_6
        return (FheInt8(bits=out_bits), overflow)

    def Mul(self, b: FheInt8) -> FheInt8:
        out = eight_bits_mul(a=self.bits, b=b.bits)
        return FheInt8(bits=out)

    def DivAndRemOverflow(self, b: FheInt8) -> (FheInt8, FheInt8, bool):
        '''
        Same as DivAndRem but also sets the overflow flag when a = -128 and b = 1. 
        Refer to `DivAndRem` for more information on signed divison. 

        Returns (quotient, remainder, div_error flag, overflow flag)

        # Overflow
        Overflow happens when a=-128 and b=-1. For more information check `DivAndRem`. 
        To check for overflow, since
            1 0 0 0 0 0 0 0 (bits of -128)
            1 1 1 1 1 1 1 1 (bit of -1)
        we negate MSB of first least signficant bits of b and AND all the bits which requires 
        at-lest 15 AND operations. 

        In case of overflow, returned quotient c = -128 and rremainder d=0
        
        '''
        (quotient, remainder, div_error) = self.DivAndRem(b)

        # overflow check
        overflow = b.bits[7] & self.bits[7]
        for i in range(7):
            overflow = overflow & (b.bits[i]&(not self.bits[i]))

        return (quotient, remainder, div_error, overflow)

    def DivAndRem(self, b: FheInt8) -> (FheInt8, FheInt8, bool):
        '''
        returns c, d, div_error s.t. 
            a = c*b + d (c: quotitent, d: remainder)
            
            AND div_error flag indicates whether division by 0 was attempted.

        If b=0, then division by 0 was attempeted and div_error flag is set to 1. In this case, 
        returned c = -1 and d = a.

        Signed division is performed by using unsigned diviaion as a subroutine. We first take 
        absolute values of a and b and then treat them as unsigned integers for divison. We then 
        adjust the sign of quotient (i.e. c) and remainder (i.e. d) as per signs of a and b: 

        |   a   |   b   |   c   |   d   |
        ---------------------------------
        |   +   |   +   |   +   |   +   |
        |   -   |   -   |   +   |   -   |
        |   +   |   -   |   -   |   +   |
        |   -   |   +   |   -   |   -   |

        Sign of `c` is obvious. For sign of `d` the key thing to note is that `d` is remainder of unsigned integer division. 
        Hence, its sign must adjusted such that a is obtained after adding d to c*b. 

        Using the table we can derive the following rules to adjust the sign: 
        (1) Change `d` to -ve whenever `a` is -ve. 
        (2) Change `c` to -ve whenever signs of `a` and `b` are different (i.e sign(a)^sign(b) == 1)

        Note that in signed representation value is -ve if MSB is set to True, otherwise +ve.

        # Caveat in quotient output when div_by_zero=1 in signed division
        To return -1 when division is attempeted by 0 we need to be careful in switching signs when 
        numerator is -ve. This is because, as seen in table, sign of quotient must be negated whenever
        numerator and denominator have opposite signs. This creates an issue when denomiator is 0 because
        0 assumes +ve sign (MSB of 0 in signed representation is 0). Hence, just negating whenever
        MSB(a)^MSB(b) == 1 also negates quotient value when division is attempeted by 0. Arbitrary unsigned division 
        returns maximum value, ie. 255, that can fit in 8 bits, whenever denomaintor = 0. Which when re-interpreted as 
        signed value is -1. If we negate -1, it changes to +1, which is incorrect output. Hence, when negating, in addition to 
        checking MSB(a)^MSB(b) == 1 we should also check div_by_zero error is False.

        # Overflow

        Signed values, unlike unsgined values, can overflow when a = -2^{N-1} and b = -1. This is because
            -2^{N-1} / -1 = 2^{N-1}
        and 2^{N-1} is out of range for signed representation using N bits. For ex, in Int8 division overflows
        whenever a = -128 and b = -1, because a/b = 128 and there's no space for 128. 

        Whenever division overflows, c = a = -2^{N-1} and d = 0. 
        '''
        pos_a = absolute(a=self.bits)
        pos_b = absolute(a=b.bits)

        div_error = is_zero(a=pos_b)

        (quotient, remainder) = arbitrary_unsigned_division(a=pos_a, b=pos_b)

        # set sign of quotient
        neg_quotient = [not i for i in quotient]
        carry = True
        for i in range(8):
            (neg_quotient[i], carry) = half_adder(A=neg_quotient[i], B=carry)
        # (self.bits[-1]^b.bits[-1]) & (not div_error) == 1 then negate quotient otherwise quotient remains unchanged
        quotient = mux_bool_vec(bit=((self.bits[-1]^b.bits[-1]) & (not div_error)), a=neg_quotient, b=quotient)

        # if (self.bits[-1]^b.bits[-1]) & (not div_error):
        #     # negate quotient
        #     quotient = [not i for i in quotient]
        #     carry = True
        #     for i in range(8):
        #         (quotient[i], carry) = half_adder(A=quotient[i], B=carry)
        # else:
        #     # quotient stays +ve
        #     pass

        # set sign of remainder
        neg_remainder = [not i for i in remainder]
        carry = True
        for i in range(8):
            (neg_remainder[i], carry) = half_adder(A=neg_remainder[i], B=carry)
        # if MSB of `a` is 1, then negate remainder otherwise remainder remains unchanged
        remainder = mux_bool_vec(bit=self.bits[-1], a=neg_remainder, b=remainder)

        # # set sign of remainder
        # if self.bits[-1]:
        #     # negate remainder if dividend is -ve
        #     remainder = [not i for i in remainder]
        #     carry = True
        #     for i in range(8):
        #         (remainder[i], carry) = half_adder(A=remainder[i], B=carry)
        # else:
        #     # remainder stays +ve
        #     pass
        
        return (FheInt8(bits=quotient), FheInt8(bits=remainder), div_error)

    def GreaterThan(self, b: FheInt8) -> bool:
        return arbitrary_signed_bit_comparator(a=self.bits, b=b.bits)

    def GreaterThanOrEqualTo(self, b: FheInt8) -> bool:
        # A>=B  = !(A<B)
        return not (self.LessThan(b))

    def LessThan(self, b: FheInt8) -> bool:
        return arbitrary_signed_bit_comparator(a=b.bits, b=self.bits)

    def LessThanOrEqualTo(self, b: FheInt8) -> bool:
        # A<=B = !(A>B)
        return not (self.GreaterThan(b=b))

    def Equals(self, b: FheInt8) -> bool:
        return arbitrary_bit_equality(a=self.bits, b=b.bits)     
class FheUint8:
    def __init__(self, bits: [bool]):
        assert len(bits) == 8
        self.bits = bits

    def from_uint8(v: int) -> FheUint8:
        # value can be in range [0, 255]
        assert v >= 0

        # extract 8 bits
        bits = []
        for i in range(0, 8):
            bits.append(((v >> i) & 1)==1)

        return FheUint8(bits=bits)

    def to_uint8(self) -> int:
        v = 0
        for i in range(8):
            v += (self.bits[i]*(1 << i))
        return v

    def Add(self, b: FheUint8) -> (FheUint8, bool):
        '''
        Adds two unsigned integers

        # Overflow
        Addition of 2 unsigned integers overflows when carry out bit, i.e. c_7, is set to 1. This is because when c_7 = 1, 
        x + y = 2^N + z = z mod{2^N}
        '''
        (out, c_7, _) = arbitrary_bit_adder(a=self.bits, b=b.bits, carry_in=False)
        return (FheUint8(bits=out), c_7)

    def Sub(self, b: FheUint8) -> (FheUint8, bool):
        '''
        Subtracts two unsigned integers

        # Overflow
        Subtraction of 2 unsigned integers overflows when carry out bit, i.e. c_7, is set to 0. This is because
            x - y = x + (2^N - y) = 2^N + x - y
        So the default case for subtraction is c_7 = 1. For overflow (or underflow, whatever you prefer) to happen, y > x. 
        If y > x, the output must be < 2^N, hence c_7 must be 0.
        '''

        (out, c_7, _) = arbitrary_bit_subtractor(a=self.bits, b=b.bits, borrow_in=False)
        return (FheUint8(bits=out), not c_7)

    def DivAndRem(self, b: FheInt8) -> (FheUint8, FheUint8, bool):
        '''
        returns c, d, div_error s.t. 
            a = c*b + d (c: quotitent, d: remainder)
            
            AND div_error flag indicates whether division by 0 was attempted.

        If b=0, then division by 0 was attempeted and div_error flag is set to 1. In this case, 
        returned c = -1 and d = a.
        '''
        div_error = is_zero(a=b.bits)
        (quotient, remainder) = arbitrary_unsigned_division(a=self.bits, b=b.bits)
        return (FheUint8(bits=quotient), FheUint8(bits=remainder), div_error)


    def Mul(self, b: FheUint8) -> FheUint8:
        out= eight_bits_mul(a=self.bits, b=b.bits)
        return FheUint8(bits=out)

    def GreaterThan(self, b: FheUint8) -> bool:
        return arbitrary_unsigned_bit_comparator(a=self.bits, b=b.bits)

    def GreaterThanOrEqualTo(self, b: FheUint8) -> bool:
        # A>=B  = !(A<B)
        return not (self.LessThan(b))

    def LessThan(self, b: FheUint8) -> bool:
        return arbitrary_unsigned_bit_comparator(a=b.bits, b=self.bits)

    def LessThanOrEqualTo(self, b: FheUint8) -> bool:
        # A<=B = !(A>B)
        return not (self.GreaterThan(b=b))

    def Equals(self, b: FheUint8) -> bool:
        return arbitrary_bit_equality(a=self.bits, b=b.bits)

def unsigned_tests():
    # Unsigned integers
    for i in range(256):
        for j in range(256):
            fhe_i = FheUint8.from_uint8(i)
            fhe_j = FheUint8.from_uint8(j)

            # Add
            want_out = (i+j)%256
            want_overflow = (i+j) != want_out
            (fhe_out, overflow) = fhe_i.Add(b=fhe_j)
            assert want_out == fhe_out.to_uint8()
            assert want_overflow == overflow

            # Sub
            want_out = (i-j)%256
            want_overflow = (i-j) != want_out
            (fhe_out, overflow) = fhe_i.Sub(b=fhe_j)
            assert want_out == fhe_out.to_uint8(), f'Want {want_out} but got {fhe_out.to_uint8()} for {i}-{j}'
            assert want_overflow == overflow, f'Want overflow {want_overflow} but got {overflow} for {i}-{j}'

            # Mul
            want_out = (i*j)%256
            fhe_out = fhe_i.Mul(b=fhe_j)
            assert want_out == fhe_out.to_uint8()

            # Division
            (fhe_quotient, fhe_remainder, div_zero_error) = fhe_i.DivAndRem(fhe_j)
            if j != 0:
                want_quotient = int(i//j)
                want_rem = i - (want_quotient*j)
                assert want_quotient == fhe_quotient.to_uint8()
                assert want_rem == fhe_remainder.to_uint8()
            else:
                assert div_zero_error

            # Comparators
            assert fhe_i.Equals(b=fhe_j) == (i==j), f'Want {i==j} but got {fhe_i.Equals(b=fhe_j)} for {i}=={j}'
            assert fhe_i.LessThan(b=fhe_j) == (i<j)
            assert fhe_i.GreaterThan(b=fhe_j) == (i>j)
            assert fhe_i.LessThanOrEqualTo(b=fhe_j) == (i<=j)
            assert fhe_i.GreaterThanOrEqualTo(b=fhe_j) == (i>=j)

def uint8_to_int8(v: int) -> int:
        assert v < 256 and v >= 0, f'v={v}'
        if v > 127:
            return -(256 - v)
        else:
            return v

def signed_tests():
    # Signed integers
    for i_unsigned in range(256):
        for j_unsigned in range(256):
            i = uint8_to_int8(i_unsigned)
            j = uint8_to_int8(j_unsigned)

            fhe_i = FheInt8.from_int8(i)
            fhe_j = FheInt8.from_int8(j)

            np_i = np.int8(i)
            np_j = np.int8(j)

            # Add
            want_out = np_i + np_j
            want_overflow = (i+j) != want_out.tolist()
            (fhe_out, overflow) = fhe_i.Add(b=fhe_j)
            assert want_out.tolist() == fhe_out.to_int8()
            assert want_overflow == overflow

            # Sub
            want_out = np_i-np_j
            want_overflow = (i-j) != want_out.tolist()
            (fhe_out, overflow) = fhe_i.Sub(b=fhe_j)
            assert want_out.tolist() == fhe_out.to_int8(), f'Want {want_out} but got {fhe_out.to_int8()} for {i}-{j}'
            assert want_overflow == overflow, f'Want overflow {want_overflow} but got {overflow} for {i}-{j}'

            # Mul
            want_out = np_i*np_j
            fhe_out = fhe_i.Mul(b=fhe_j)
            assert want_out.tolist() == fhe_out.to_int8()

            
            if i == -128 and j == -1:
                # Division overflow
                (fhe_quotient, fhe_remainder, div_zero_error, overflow) = fhe_i.DivAndRemOverflow(fhe_j)
                assert fhe_quotient.to_int8() == -128, f'Want -128 but got {fhe_quotient.to_int8()}'
                assert fhe_remainder.to_int8() == 0
                assert div_zero_error == False
                assert overflow == True

                (fhe_quotient, fhe_remainder, div_zero_error) = fhe_i.DivAndRem(fhe_j)
                assert fhe_quotient.to_int8() == -128
                assert fhe_remainder.to_int8() == 0
                assert div_zero_error == False
            else:
                (fhe_quotient, fhe_remainder, div_zero_error) = fhe_i.DivAndRem(fhe_j)
                if j != 0:
                    want_quotient = int(np_i/np_j)
                    want_rem = i - (want_quotient*j)
                    assert want_quotient == fhe_quotient.to_int8(), f'Want {want_quotient} but got {fhe_quotient.to_int8()} for {i}/{j}'
                    assert want_rem == fhe_remainder.to_int8()
                    assert div_zero_error == False
                else:
                    assert div_zero_error
                    assert fhe_quotient.to_int8() == -1
                    assert fhe_remainder.to_int8() == i

            # Comparators
            assert fhe_i.Equals(b=fhe_j) == (i==j), f'Want {i==j} but got {fhe_i.Equals(b=fhe_j)} for {i}=={j}'
            assert fhe_i.LessThan(b=fhe_j) == (i<j)
            assert fhe_i.GreaterThan(b=fhe_j) == (i>j)
            assert fhe_i.LessThanOrEqualTo(b=fhe_j) == (i<=j)
            assert fhe_i.GreaterThanOrEqualTo(b=fhe_j) == (i>=j)

# a = FheInt8.from_int8(-128)
# b = FheInt8.from_int8(-1) 
# # c = a.Mul(b)
# # print(c.to_int8())
# (c, c_rem, div_error, overflow) = a.DivAndRemOverflow(b)
# print(c.to_int8(), c_rem.to_int8(), div_error, overflow)
# a = FheUint8.from_uint8(10)
# b = FheUint8.from_uint8(11)
# print(a.GreaterThan(b))
# print(a.LessThan(b))
signed_tests()