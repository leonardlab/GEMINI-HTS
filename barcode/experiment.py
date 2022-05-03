
from .classifier import BarcodeClassifier
from .descriptor import BarcodeDescriptor

from segments.experiment import Experiment
from segments.util import reverse_complement


class BarcodeExperiment(Experiment):

    def aliases(self):
        return [ "barcode" ]

    def fragmentDescriptor(self):
        return BarcodeDescriptor().fragmentDescriptor()

    def classifier(self, descriptor):
        return BarcodeClassifier(descriptor)

    @property
    def knownBarcodes(self):
        return [
            "TATAACTACCTGCGGTGACG",
            "TCACCTGGCTAGACTCAACC",
            "ATCAATGCGGACACGATGTG",
            "GGATGCCCAAGTTCTGAGTA",
        ]
