
from segments.descriptor import Descriptor, M


# 2022_07_30: sequences added to both fp/rp
# 2022_08_29: fp and rp extended by 4 nt on each side, looks like fp has a toggle change

class BarcodeDescriptor(Descriptor):

    def __init__(self):
        Descriptor.__init__(self)
        self.maxAllowedErrors = 4
        self.handleIndels = False
        self.barcodeLen = 20
        #fp = "TTTTTTAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        fp  = "TTTTTTAATGATACGGCGACCACCGAGATATACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        self.fp = fp + "TCCGTCACGCCAATCATAGG"
        self.fp = "TTTT" + self.fp
        #rp = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTTCTATCTCGTATGCCGTCTTCTGCTTG"
        rp = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
        rp = rp + "AAAA"
        self.rp = "GAAACGGTATTCGCACACGAGA" + rp
        self.left = "AACTGCACGCGT"
        self.right = "GAGCTCTACTGT"

    def fragmentDescriptor(self):
        fp = self.makeSegd("fPrimer", self.fp, index = True, fixedLength = True,
                           maxAllowedErrors = self.maxAllowedErrors, handleIndels = self.handleIndels)
        l = self.makeSegd("left", self.left, index = True, fixedLength = True,
                          maxAllowedErrors = self.maxAllowedErrors, handleIndels = self.handleIndels)
        b = self.makeSegd("barcode", wildcard = True, fixedLength = self.barcodeLen)
        r = self.makeSegd("right", self.right, index = True, fixedLength = True,
                          maxAllowedErrors = self.maxAllowedErrors, handleIndels = self.handleIndels)
        rp = self.makeSegd("rPrimer", self.rp, index = True, fixedLength = True,
                           maxAllowedErrors = self.maxAllowedErrors, handleIndels = self.handleIndels)
        return M.SequenceDescriptor(segments = [ fp, l, b, r, rp ])
