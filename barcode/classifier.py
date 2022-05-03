
from segments.classifier import Classifier, Failures
from segments.util import Indel, fs, writeJsonToPath, reverse_complement

from ngs import model
BarcodeResult = model.NGSResult


class BarcodeClassifier(Classifier):

    def __init__(self, descriptor):
        Classifier.__init__(self, descriptor)
        self.minimumRequiredOverlap = 10

    def count(self, fr):
        fragment, sres, res, r1s, r2s = fr.f, fr.s, fr.c, fr.r1, fr.r2
        assert(res and res.category)

        ctrs = self.counters
        ctrs.countedPairs += 1
        ctrs.incrementKey("category:{}".format(res.category))

        if fragment.success:
            ctrs.incrementKey("overlap:{}".format(fragment.overlap))

        if res.success:
            ctrs.matchedPairs += 1
        else:
            if not fragment.success:
                ctrs.incrementKey("fragmentFailure:{}".format(fragment.failure))
            elif sres and sres.failure:
                ctrs.incrementKey("segmentFailure:{}".format(sres.failure))
                if sres.failedSegment:
                    ctrs.incrementKey("failedSegment:{}:{}".format(sres.failedSegment, sres.failure))
            return

        ctrs.incrementKey("errors:{}".format(res.errorCount))

        barcode = sres.subsequenceForSegment(sres.segmentWithKey('barcode'))
        if res.errorCount:
            ctrs.matchWithErrors += 1
            ctrs.incrementKey("barcodeWithErrors:{}".format(barcode))
        else:
            ctrs.perfectMatch += 1
            ctrs.incrementKey("barcode:{}".format(barcode))

        fp, rp = sres.segmentWithKey("fPrimer"), sres.segmentWithKey("rPrimer")
        if fp and fp.errorCount():
            ctrs.fpError += 1
            ctrs.incrementKey("fpErr:{}".format(fp.errorCount()))
        if rp and rp.errorCount():
            ctrs.rpError += 1
            ctrs.incrementKey("rpErr:{}".format(rp.errorCount()))

    def _classifyBasic(self, fragment, sres):
        res = BarcodeResult(success = fragment.success and sres.success)
        if res.success and fragment.overlap <= self.minimumRequiredOverlap:
            res.category = "overlapFailure"
            res.success = False
            return res
        if res.success:
            res.category = None
        else:
            if not fragment.success:
                assert(fragment.failure)
                res.category = fragment.failure
            elif not sres:
                res.category = "segmentFailure"
            elif sres.failure:
                res.category = sres.failure
            else:
                ASSERT_NOT_REACHED
            assert(res.category)
        return res

    def classify(self, fragment, sres):
        res = self._classifyBasic(fragment, sres)
        if res.success:
            self._classifyPrimers(fragment, sres, res)
            res.category = res.category or "fullMatch"
        assert(res.category)
        return res

    def _classifyPrimers(self, fragment, sres, res):
        fp, rp = sres.segmentWithKey("fPrimer"), sres.segmentWithKey("rPrimer")
        if fp and fp.errorCount():
            res.category = res.category or "primerError"
        if rp and rp.errorCount():
            res.category = res.category or "primerError"
