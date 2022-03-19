
from segments.classifier import Classifier, Failures
from segments.util import Indel, fs, writeJsonToPath, reverse_complement
from . import model as M


LARGE_INDEL_SIZE = 4

class NGSClassifier(Classifier):

    def __init__(self, descriptor):
        Classifier.__init__(self, descriptor)
        self.skipBarcodePruning = False
        self.barcodeTracking = False
        self.minimumRequiredOverlap = 10
        self.reads = []

    def _addRead(self, res, sequenceId):
        read = res.read
        assert(read)
        read.mutations = []
        for m in (res.mutations or []):
            read.mutations.append(M.NGSMutation(mutType = fs.mutation, site = m))
        for i in (res.indels or []):
            m = M.NGSMutation(site = i.errorIndex, size = i.size)
            if i.substSize:
                m.mutType = fs.substitution
                m.substSize = i.substSize
            elif i.insertType:
                m.mutType = fs.insertion
            else:
                m.mutType = fs.deletion
            read.mutations.append(m)
        assert(read.mutations)
        self.reads.append(read)

    def saveReads(self, key):
        fn = "{}_reads.json".format(key)
        writeJsonToPath({ "reads" : [ x.json(skipTypes = True) for x in self.reads ] }, fn)
        self.info("Saved reads to:", fn)

    def getBarcode(self, fr):
        return  fr.r1.characters[:5] + "_" + reverse_complement(fr.r2.characters[:5])

    def count(self, fr):
        fragment, sres, res, r1s, r2s = fr.f, fr.s, fr.c, fr.r1, fr.r2
        assert(res and res.category)

        ctrs = self.counters
        ctrs.countedPairs += 1
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

        ctrs.incrementKey("category:{}".format(res.category))
        ctrs.incrementKey("errors:{}".format(res.errorCount))

        if res.errorCount:
            self._addRead(res, r1s.sequenceId)
            ctrs.matchWithErrors += 1
        else:
            ctrs.perfectMatch += 1

        # TODO: DELME, since we have reads now? or nbd? well, showCounts is nice.
        muts = res.mutations or []
        indels = res.indels or []

        for mut in muts:
            ctrs.incrementKey("err:{}".format(mut))
            ctrs.incrementKey("mut:{}".format(mut))

        for i in indels:
            ctrs.incrementKey("err:{}".format(i.errorIndex))
            if i.substSize:
                ctrs.incrementKey("subst:{}".format(i.errorIndex))
                ctrs.incrementKey("substSize:{}_{}".format(i.size, i.substSize))
                ctrs.incrementKey("substAtSize:{}:{}_{}".format(i.errorIndex, i.size, i.substSize))
                for ei in range(i.substSize):
                    ctrs.incrementKey("goner:{}".format(i.errorIndex + ei))
                if (i.substSize > LARGE_INDEL_SIZE) or (i.size > LARGE_INDEL_SIZE):
                    ctrs.incrementKey("largeSubst")
            elif i.insert_type:
                ctrs.incrementKey("ins:{}".format(i.errorIndex))
                ctrs.incrementKey("insSize:{}".format(i.size))
                ctrs.incrementKey("insAtSize:{}:{}".format(i.errorIndex, i.size))
                if i.size >= LARGE_INDEL_SIZE:
                    ctrs.incrementKey("largeIns")
            else:
                ctrs.incrementKey("del:{}".format(i.errorIndex))
                ctrs.incrementKey("delSize:{}".format(i.size))
                ctrs.incrementKey("delAtSize:{}:{}".format(i.errorIndex, i.size))
                for ei in range(i.size):
                    ctrs.incrementKey("goner:{}".format(i.errorIndex + ei))
                if i.size >= LARGE_INDEL_SIZE:
                    ctrs.incrementKey("largeDel")

        if muts:
            ctrs.incrementKey("muts:{}".format(len(muts)))

        fp, rp = sres.segmentWithKey("fPrimer"), sres.segmentWithKey("rPrimer")
        if fp and fp.errorCount():
            ctrs.fpError += 1
            ctrs.incrementKey("fpErr:{}".format(fp.errorCount()))
        if rp and rp.errorCount():
            ctrs.rpError += 1
            ctrs.incrementKey("rpErr:{}".format(rp.errorCount()))

    def _classifyBasic(self, fragment, sres):
        res = M.NGSResult(success = fragment.success and sres.success)
        res.read = M.NGSRead(sequenceId = fragment.seq.sequenceId, success = res.success, overlap = fragment.overlap)
        if res.success and fragment.overlap <= self.minimumRequiredOverlap:
            res.category = "overlapFailure"
            res.success = False
            return res
        if res.success:
            res.category = None
        else:
            if not fragment.success:
                assert(fragment.failure)
                res.category = "fragmentFailure"
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
            muts, indels = self._targetErrors(fragment, sres)
            self._classifyErrors(muts, indels, res)
            self._classifyPrimers(fragment, sres, res)
            ts = self._targetSegment(sres)
            res.targetDescription = self._makeTargetDescription(ts)
            res.targetSequence = sres.subsequenceForSegment(ts)
            res.category = res.category or "fullMatch"
        assert(res.category)
        return res

    def _makeTargetDescription(self, sr):
        desc = []
        if sr.errors:
            desc.append("![{}]".format(",".join([str(x) for x in sr.errors])))
        for ind in sr.indels:
            desc.append("{}".format(ind))
        return "".join(desc)

    def _targetSegment(self, sres):
        return sres.segmentWithKey("target")

    def _targetErrors(self, fragment, sres):
        t = sres.segmentWithKey("target")
        assert(t)
        muts = sorted(t.errors or [])
        indels = sorted(t.indels or [], key = lambda i : i.errorIndex)
        return muts, indels

    def _classifyErrors(self, muts, indels, res):
        desc = "m"
        for mut in muts:
            desc += "!{}".format(mut)
        for i in indels:
            if i.substSize:
                desc += "+{}-{}@{}".format(i.size, i.substSize, i.errorIndex)
            elif i.insert_type:
                desc += "+{}@{}".format(i.size, i.errorIndex)
            else:
                desc += "-{}@{}".format(i.size, i.errorIndex)
        res.description = "full" if desc == "m" else desc
        res.mutations = muts or None
        res.mutationCount = len(muts)
        res.indels = indels or None
        res.errorCount = res.mutationCount + len(indels or [])

    def _classifyPrimers(self, fragment, sres, res):
        fp, rp = sres.segmentWithKey("fPrimer"), sres.segmentWithKey("rPrimer")
        if fp and fp.errorCount():
            res.category = res.category or "primerError"
        if rp and rp.errorCount():
            res.category = res.category or "primerError"
