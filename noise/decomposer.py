from __future__ import annotations

class Decomposer():
    def __init__(self, d_a: Integer, d_b: Integer, logQ: Integer, logB: Integer, max_logB: Integer = None):
        self.d_a = d_a
        self.d_b = d_b
        self.logB = logB
        if max_logB == None:
            max_logB = logQ
        self.max_logB = max_logB
        self.logQ = logQ

    def ignore_bits_a(self):
        return max(0, self.logQ - (self.d_a*self.logB))

    def ignore_bits_b(self):
        return max(0, self.logQ - (self.d_b*self.logB))

    def optimize(self, f, log_fail: float) -> float:
        pr_fail = f()

        while pr_fail < log_fail and self.logB < self.max_logB:
            self.logB += 1
            pr_fail = f()
            if pr_fail > log_fail:
                self.logB -= 1
                pr_fail = f()
                break

        while pr_fail <= log_fail and self.d_a != 1:
            self.d_a -= 1
            pr_fail = f()
            if pr_fail > log_fail:
                self.d_a += 1
                pr_fail = f()
                break

        if self.d_b != None:
            while pr_fail <= log_fail and self.d_b != 1:
                self.d_b -= 1
                pr_fail = f()
                if pr_fail > log_fail:
                    self.d_b += 1
                    pr_fail = f()
                    break

        return pr_fail
    
    @staticmethod
    def double_decomposer(d_a: Integer, d_b: Integer, logQ: Integer, logB: Integer, max_logB: Integer = None):
        return Decomposer(d_a=d_a, d_b=d_b, logQ=logQ, logB=logB, max_logB=max_logB)
    
    @staticmethod
    def single_decomposer(d: Integer, logQ: Integer, logB: Integer, max_logB: Integer = None):
        return Decomposer(d_a=d, d_b=None, logB=logB, logQ=logQ, max_logB=max_logB)