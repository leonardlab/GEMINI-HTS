
from segments.descriptor import Descriptor, M


class BarcodeDescriptor(Descriptor):

    def __init__(self):
        Descriptor.__init__(self)
        self.maxAllowedErrors = 4
        self.handleIndels = False
        self.barcodeLen = 20
        self.fp = "TCCGTCACGCCAATCATAGG"
        self.rp = "GAAACGGTATTCGCACACGAGA"
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
