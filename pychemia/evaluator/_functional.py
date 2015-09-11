import os
import time
import numpy as np
from threading import Thread


class ObjectiveFunction:
    def __init__(self):
        self.population = None

    def initialize(self, population):
        self.population = population

    def ids_sorted(self, selection):
        values = np.array([self.population.value(i) for i in selection])
        argsort = np.argsort(values)
        return np.array(selection)[argsort]

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.population.value(i)
        return ret


class Evaluator:
    def __init__(self):
        self.process = None
        self.thread = None
        self.population = None

    @property
    def is_running(self):
        if self.thread is not None:
            return self.thread.is_alive()
        else:
            return False

    def initialize(self, population):
        self.population = population

    def evaluate(self, i):
        x = self.population.coordinate(i)
        y = self.population.function(x)
        if y is not None:
            self.population.set_value(i, y)
            self.population.evaluated.append(i)

    def run(self):

        def worker(evaluator, population):
            while True:
                for i in population.actives:
                    if not population.is_evaluated(i):
                        evaluator.evaluate(i)
                time.sleep(1)
                if os.path.exists('stop'):
                    os.remove('stop')
                    return

        # self.process = Process(target=worker, args=(self.population,))
        # self.process.start()
        self.thread = Thread(target=worker, args=(self, self.population,))
        self.thread.daemon = True
        self.thread.start()

    def stop(self):
        if self.process is not None and not self.process.is_alive():
            self.process.terminate()
