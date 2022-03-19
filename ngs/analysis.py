
from segments.util import fs, jsonAtPath, writeCsvToPath

from . import model as M
from .experiment import NGSExperiment
from .run import NGSRunner



class NGSAnalysis(object):

    def __init__(self):
        self.experiment = NGSExperiment()
        self.keys = self.ngs_experiment2_keys()

    def ngs_experiment2_keys(self, short65k = False):
        res = [ "s{}_{}".format(i + 1, "65k" if int(short65k) else "full") for i in range(5) ]
        if short65k:
            res = [ res[0], res[2], res[4] ]
        return res

    def run(self):
        NGSRunner(self.experiment, self.keys).runAll()
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

    def parseReads(self, key):
        reads = jsonAtPath("{}_reads.json".format(key))
        counts = { }
        for r in reads['reads']:
            overlaps = counts.setdefault('overlap', {})
            o = r.get('overlap', 0)
            overlaps[o] = overlaps.get(o, 0) + 1
            for mi in r['mutations']:
                m = M.NGSMutation().loadFromJson(mi)
                sites = counts.setdefault(m.mutType, {})
                delSize = 0
                if m.mutType == fs.deletion:
                    delSize = m.size
                elif m.mutType == fs.substitution:
                    delSize = m.substSize
                site = m.site + (delSize >> 1)
                sites[site] = sites.get(site, 0) + 1
                errSites = counts.setdefault(fs.error, {})
                errSites[site] = errSites.get(site, 0) + 1
                if not delSize:
                    continue
                gonerSites = counts.setdefault(fs.goner, {})
                for i in range(delSize):
                    gonerSites[m.site + i] = gonerSites.get(m.site + i, 0) + 1
        return counts

    def analysis(self):
        self.errorsFromReads()
        self.results()
        self.insertSizes()
        self.delSizes()
        self.scatter('insertion')
        self.scatter('deletion')
        self.scatter('substitution')

    def errorsFromReads(self):
        keys = self.keys
        rows = [ [ "", "" ] + keys ]
        allCounts = { }
        cnts = [ self.counts(k) for k in keys ]
        x = self.experiment
        tgt = x.target
        includePctOfMatched = False
        includeRatioToFirst = False

        for k in keys:
            allCounts[k] = self.parseReads(k)

        for etype in [ fs.error, fs.mutation, fs.insertion, fs.deletion, fs.substitution, fs.goner ]:
            estr = fs.toString(etype)
            rows.append([f"{estr}:"] + ([""] * (3 * len(keys) + 3)))
            for idx in range(len(tgt)):
                row = [ idx, tgt[idx] ] + [ allCounts[k].get(etype, {}).get(idx, 0) for k in keys ]
                if includePctOfMatched:
                    row += [ "" ] + [ "{:.2f}%".format(100 * allCounts[k].get(etype, {}).get(idx, 0) / cnt['matchedPairs']) for k, cnt in zip(keys, cnts) ]
                if includeRatioToFirst:
                    row += [ "" ] + [ "{:.4f}".format(allCounts[k].get(etype, {}).get(idx, 0) / allCounts[keys[0]].get(etype, {}).get(idx, 0)) if allCounts[keys[0]].get(etype, {}).get(idx, 0) else "--" for k in keys ]
                rows.append(row)

        rows.append([""] * (len(keys) + 1))
        for o in range(100):
            rows.append([ f"overlap_{o}", ' ' ] + [ allCounts[k].get('overlap', {}).get(o, 0) for k in keys ])

        fn = "read_sample_errors.csv"
        writeCsvToPath(rows, fn)
        print(f"Wrote {fn}")

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
        _do('matchedPairs')
        uniques = _do('uniqueBarcode', pcts = [ 1 << 20 ] * 5, pctTitle = '% of max possible')
        _do('barcodeSkipped', pcts = totals)
        #_do('barcodeMismatch', pcts = uniques, pctTitle = '% of uniques')
        _do('topBarcode')

        rows.append([ "" for x in rows[0] ])

        counted = _do('countedPairs', pcts = totals)
        #_do('uniqueFragment')
        _do('matchWithErrors', pcts = counted)
        errs = [ c.get('errors', {}) for c in splits ]
        rows.append([ 'avgErrors' ] + [ "{:.4f}".format(sum([ int(k) * int(v) for k, v in e.items()]) / max(1, sum([v for k, v in e.items() if int(k)]))) for e in errs ])
        rows.append([ 'avgErrorsInclZero' ] + [ "{:.4f}".format(sum([ int(k) * int(v) for k, v in e.items()]) / t) for e, t in zip(errs, counted) ])

        rows.append([ "" for x in rows[0] ])

        _do('largeDel')
        _do('largeIns')
        _do('largeSubst')
        _do('fpError')
        _do('rpError')
        _do('fragmentReverseOrder')

        rows.append([ "" for x in rows[0] ])

        cats = set()
        for s in splits:
            for c in s.get('category', {}).keys():
                cats.add(c)
        for cat in cats:
            _do(cat, pcts = counted, vals = [s.get('category', {}) for s in splits])

        failures = set()
        for s in splits:
            ff = { k : v for k, v in s.get('fragmentFailure', {}).items() }
            sf = { k : v for k, v in s.get('segmentFailure', {}).items() }
            ff.update(sf)
            s['failures'] = ff
            for f in ff.keys():
                failures.add(f)
        for f in failures:
            _do(f, pcts = counted, vals = [s.get('failures', {}) for s in splits])

        writeCsvToPath(rows, "sample_counts.csv")
        print("Wrote results to sample_counts.csv")

    def insertSizes(self):
        keys = self.keys
        rows = [ [ "" ] + keys ]
        allCounts = { }
        maxSize = 0
        for k in keys:
            reads = jsonAtPath("{}_reads.json".format(k))
            counts = {}
            for r in reads['reads']:
                for m in r['mutations']:
                    if m['mutType'] == fs.insertion or m['mutType'] == fs.substitution:
                        counts[m['size']] = counts.get(m['size'], 0) + 1
                        maxSize = max(maxSize, m['size'])
            allCounts[k] = counts

        for sz in range(1, maxSize + 1):
            rows.append([ str(sz) ] + [ allCounts[k].get(sz, 0) for k in keys ])
        fn = "insert_sizes.csv"
        writeCsvToPath(rows, fn)
        print(f"Wrote {fn}")

    def delSizes(self):
        keys = self.keys
        rows = [ [ "" ] + keys ]
        allCounts = { }
        maxSize = 0
        for k in keys:
            reads = jsonAtPath("{}_reads.json".format(k))
            counts = {}
            for r in reads['reads']:
                for m in r['mutations']:
                    sz = 0
                    if m['mutType'] == fs.deletion:
                        sz = m['size']
                    elif m['mutType'] == fs.substitution:
                        sz = m['substSize']
                    counts[sz] = counts.get(sz, 0) + 1
                    maxSize = max(maxSize, sz)
            allCounts[k] = counts

        for sz in range(1, maxSize + 1):
            rows.append([ str(sz) ] + [ allCounts[k].get(sz, 0) for k in keys ])
        fn = "del_sizes.csv"
        writeCsvToPath(rows, fn)
        print(f"Wrote {fn}")

    def scatter(self, errorType = None):
        for key in self.keys:
            self.sizeScatter(key, errorType)

    def sizeScatter(self, key, errorType = None):
        try:
            from .plot import scatter
        except:
            print("Please install matplotlib to generate scatterplots")
            return

        errorType = getattr(fs, errorType)
        reads = jsonAtPath("{}_reads.json".format(key))

        counts = {}
        def addSiteSize(st, sz):
            c = counts.setdefault(st, {})
            c[sz] = c.get(sz, 0) + 1

        for r in reads['reads']:
            for mi in r['mutations']:
                m = M.NGSMutation().loadFromJson(mi)
                if errorType == fs.deletion:
                    site = m.site + (m.size >> 1)
                    if m.mutType == fs.deletion:
                        addSiteSize(m.site + (m.size >> 1), m.size)
                    #elif m.mutType == fs.substitution:
                    #    addSiteSize(m.site + (m.substSize >> 1), m.substSize)
                elif errorType == fs.substitution:
                    if m.mutType == fs.substitution:
                        addSiteSize(m.site + (m.substSize >> 1), 0 - m.substSize)
                        addSiteSize(m.site + (m.substSize >> 1), m.size)
                else:
                    assert(errorType == fs.insertion)
                    if m.mutType == fs.insertion:
                        addSiteSize(m.site, m.size)
                    #elif m.mutType == fs.substitution:
                    #    addSiteSize(m.site + (m.substSize >> 1), m.size)

        import math
        xs, ys, zs, cols = [], [], [], [] if fs.deletion == errorType else None
        bases = [ 26, 27, 28 ]
        for site, szs in counts.items():
            for sz, count in szs.items():
                xs.append(site)
                ys.append(sz)
                zs.append(4 * math.sqrt(count))
                if fs.deletion == errorType:
                    base = site - (sz >> 1)
                    if (base in bases) or (base + sz in bases):
                        cols.append('red')
                    else:
                        cols.append('blue')

        ymin = -30 if fs.substitution == errorType else 0
        fig = scatter(xs, ys, sizes = zs, show = False, xlim = [0, 64], ylim = [ymin, 50], cols = cols)
        fn = "{}_scatter_{}".format(fs.toString(errorType), key)
        for ext in [ 'png', 'pdf' ]:
            saveName = "{}.{}".format(fn, ext)
            fig.savefig(saveName)
            print("Saved {}".format(saveName))

        if not xs:
            print("Nothing to scatter: {}".format(key))
            return
        fn = "{}_scatter_{}.csv".format(fs.toString(errorType), key)
        maxSite = max(64, max(xs) + 1)
        maxSize = max(ys) + 1
        minSize = min(0, min(ys))
        rows = [ [ "" ] + [ "size_{}".format(x) for x in range(minSize, maxSize) ] ]
        for site in range(maxSite):
            rows.append([ "site_{}".format(site) ] + [ str(counts.get(site, {}).get(sz, 0)) for sz in range(minSize, maxSize) ])
        writeCsvToPath(rows, fn)
        print("  (data in {})".format(fn))


def pctfmt(x, base):
    return "{:.2f}%".format(100.*x/base) if base else "--"
