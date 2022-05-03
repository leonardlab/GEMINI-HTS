import os
import subprocess

from segments.util import fs, jsonAtPath, movePath, writeCsvToPath
from segments.runner import Runner

from .experiment import BarcodeExperiment


class BarcodeAnalysis(object):

    def __init__(self, debug = False):
        self._debug = debug
        self.experiment = BarcodeExperiment()
        self.keys = self.experimentKeys()

    def experimentKeys(self):
        short65k = self._debug
        self._samples = 8
        res = [ "s{}_{}".format(i + 1, "65k" if int(short65k) else "full") for i in range(self._samples) ]
        if short65k:
            res = [ res[0], res[2], res[4] ]
        return res

    def prepFastq(self):
        if os.path.exists('s{}_full_R2.fastq'.format(self._samples)):
            return
        for item in os.listdir('.'):
            if not item.endswith('.fastq.gz'):
                continue
            if not item.startswith('00'):
                continue
            sample = int(item[:3])
            assert(sample >=0  and  sample <= 8)
            if '_R1_' in item:
                movePath(item, 's{}_full_R1.fastq.gz'.format(sample))
            else:
                assert('_R2_' in item)
                movePath(item, 's{}_full_R2.fastq.gz'.format(sample))
        for sample in range(1, self._samples + 1):
            for r in range(1, 3):
                f = 's{}_full_R{}.fastq'.format(sample, r)
                if not os.path.exists(f):
                    assert(os.path.exists(f + '.gz'))
                    print("Decompress: {}".format(f))
                    subprocess.check_call([ 'gzip', '-d', f + '.gz' ])
                assert(os.path.exists(f))
            short = 's{}_65k_R2.fastq'.format(sample)
            if not os.path.exists(short):
                assert(os.path.exists("shortify.sh")) # copy it in
                subprocess.check_call([ "./shortify.sh", f"s{sample}" ])
            assert(os.path.exists(short))

    def run(self):
        if self._debug:
            return self.debug()
        self.prepFastq()
        Runner(self.experiment, self.keys).runAll()
        self.analysis()

    def counts(self, key):
        return jsonAtPath("{}_counts.json".format(key))

    def splitCounts(self, c):
        s = {}
        for k, v in c.items():
            cur = s
            bits = k.split(':')
            subs, last = bits[:-1], bits[-1]
            if not subs:
                last = "_{}".format(last) # b/c some things are both lone keys and subkeys
            for sk in subs:
                subk = sk
                if subk not in cur:
                    cur[subk] = {}
                cur = cur[subk]
            cur[last] = v
        return s

    def analysis(self):
        self.results()

    def results(self):
        keys = self.keys
        rows = [ [""] + keys ]
        counts = [ self.counts(k) for k in self.keys ]
        splits = [ self.splitCounts(c) for c in counts ]


        def _do(key, desc = None, pcts = None, pctTitle = None, vals = None):
            vals = [ c.get(key, 0) for c in vals or counts ]
            rows.append([ desc or key ] + vals)
            if pcts:
                rows.append([ pctTitle or '' ] + [ pctfmt(r, x) for r, x in zip(vals, pcts) ])
            return vals

        totals = _do('totalPairs')
        matched = _do('matchedPairs', pcts = totals)
        _do('matchWithErrors', pcts = totals)
        rows.append([ "" for x in rows[0] ])

        rows.append([ "barcodes: " ] + [ "" for x in rows[0][1:] ])
        bcs = [ s.get('barcode', {}) for s in splits ]
        known = self.experiment.knownBarcodes
        for kb in known:
            _do(kb, pcts = matched, vals = bcs)
        lens = [ { "unknown" : len([bck for bck in bc.keys() if bck not in known]) } for bc in bcs ]
        _do("unknown", pcts = matched, vals = lens)
        rows.append([ "" for x in rows[0] ])

        rows.append([ "categories: " ] + [ "" for x in rows[0][1:] ])
        cats = set()
        for s in splits:
            for c in s.get('category', {}).keys():
                cats.add(c)
        cats = sorted(cats, reverse = True, key = lambda c : sum([s.get('category', {}).get(c, 0) for s in splits ]))
        for cat in cats:
            _do(cat, pcts = totals, vals = [s.get('category', {}) for s in splits])

        rows.append([ "" for x in rows[0] ])

        writeCsvToPath(rows, "sample_counts.csv")
        print("Wrote results to sample_counts.csv")


    def debug(self):
        Runner(self.experiment, self.keys).runSingle(diagram = True, allDiagrams = False)


def pctfmt(x, base):
    return "{:.2f}%".format(100.*x/base) if base else "--"

def main():
    ba = BarcodeAnalysis(debug = False)
    ba.run()
