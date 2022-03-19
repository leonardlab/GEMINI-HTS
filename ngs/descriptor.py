
from segments.descriptor import Descriptor, M


class NGSDescriptor(Descriptor):

    def __init__(self):
        Descriptor.__init__(self)
        self.indexLen = 5
        self.maxAllowedPrimerErrors = 0 # LL EXP 6: 4
        self.fp = "GAGAAGCATGACGGACAAGTACAG"
        self.rp = "GTGGCAAACTGGTACTTTGGGA"
        self.target = "GCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCCTTCTGGGCAGTTGATGCC"

    def fragmentDescriptor(self):
        li = self.makeSegd("lIndex", wildcard = True, fixedLength = 5)
        fp = self.makeSegd("fPrimer", self.fp, index = True, fixedLength = True,
                           maxAllowedErrors = self.maxAllowedPrimerErrors, handleIndels = True)
        t = self.makeSegd("target", self.target, index = True, fixedLength = True,
                          maxAllowedErrors = 4, handleIndels = True)
        rp = self.makeSegd("rPrimer", self.rp, index = True, fixedLength = True,
                           maxAllowedErrors = self.maxAllowedPrimerErrors, handleIndels = True)
        ri = self.makeSegd("rIndex", wildcard = True, fixedLength = 5)
        return M.SequenceDescriptor(segments = [ li, fp, t, rp, ri ])
