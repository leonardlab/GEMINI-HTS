import multiprocessing
import time

from segments.fragment_set import FragmentSet
from segments.processor import FragmentSetProcessor


class NGSRunner(object):

    def __init__(self, exp, keys):
        self.experiment = exp
        self.keys = keys

    def _run(self, key):
        proc = FragmentSetProcessor(self.experiment, FragmentSet(key))
        proc.barcodePrune = False
        proc.autosave()
        proc.run()
        proc.saveResults()
        proc.classifier.saveReads(key)

    def _workerDispatch(self, maxWorkers = 8):
        from segments.util import _speedup
        assert(_speedup)
        workers = []
        for key in self.keys:
            worker = multiprocessing.Process(target = self._run, args = (key,))
            workers.append(worker)
            worker.start()
            while len(workers) >= maxWorkers:
                time.sleep(1)
                workers = [ w for w in workers if w.is_alive() ]
        for worker in workers:
            worker.join()

    def runAll(self):
        self._workerDispatch()
