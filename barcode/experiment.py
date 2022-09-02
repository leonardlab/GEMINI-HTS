
from .classifier import BarcodeClassifier
from .descriptor import BarcodeDescriptor

from segments.experiment import Experiment
from segments.util import reverse_complement

# 2022_07_30: some new barcodes added
# 2022_08_29: 6 different barcodes supplied

class BarcodeExperiment(Experiment):

    def aliases(self):
        return [ "barcode" ]

    def fragmentDescriptor(self):
        return BarcodeDescriptor().fragmentDescriptor()

    def classifier(self, descriptor):
        return BarcodeClassifier(descriptor)

    @property
    def knownBarcodes(self):
        if 1:
            return [
                "TATAACTACCTGCGGTGACG",
                "TCACCTGGCTAGACTCAACC",
                "ATCAATGCGGACACGATGTG",
                "GGATGCCCAAGTTCTGAGTA",
                "GCGCTAATCGGTCCTCCCGA",
                "CCGACATCCAACTGATCAAG",
            ]
        return [
            "TATAACTACCTGCGGTGACG",
            "TCACCTGGCTAGACTCAACC",
            "ATCAATGCGGACACGATGTG",
            "GGATGCCCAAGTTCTGAGTA",
            "GCGCTAATCGGTCCTCCCGA",
            "CCGACATCCAACTGATCAAG",
            "GGCCCGAATATAATTCGGAT",
            "AGCTCGTTACTCAGCACCTA",
            "GTGAGGGCTCTACGTACGAG",
            "TCGCAGTCCATGGGCAGACG",
        ]
