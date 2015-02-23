__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta, abstractmethod
from pychemia import log
import time


class Searcher():
    __metaclass__ = ABCMeta
    """
    Abstract class for all optimization algorithms that uses fixed and blocked
    Generations
    """

    def __init__(self, population, fraction_evaluated=1.0, generation_size=32, stabilization_limit=10):
        self.population = population
        self.generation_size = generation_size
        self.fraction_evaluated = fraction_evaluated
        self.stabilization_limit = stabilization_limit
        self.current_generation = 0
        self.generation = {}

    def enforce_generation_size(self):

        while True:
            size_population = len(self.population.members)
            size_active = len(self.population.actives)
            size_evaluated = len(self.population.evaluated)
            size_actives_evaluated = len(self.population.actives_evaluated)

            if size_active == self.generation_size:
                break
            elif size_active > self.generation_size:
                # Overpopulated removing some members
                if size_actives_evaluated >= self.generation_size:
                    candidates = self.population.ids_sorted(self.population.actives_evaluated)
                    index = 0
                    while len(self.population.actives) > self.generation_size:
                        if self.population.actives[index] not in candidates:
                            self.population.disable(self.population.actives[index])
                        index += 1
                else:
                    for i in range(self.generation_size - size_actives_evaluated):
                        self.population.disable(self.population.actives_no_evaluated[i])
            else:
                if size_evaluated >= self.generation_size:
                    candidates = self.population.ids_sorted(self.population.evaluated)
                    index = 0
                    while len(self.population.actives) < self.generation_size:
                        if self.population.evaluated[index] not in self.population.actives:
                            self.population.enable(self.population.evaluated[index])
                        index += 1
                else:
                    self.population.random_population(self.generation_size - size_active)

    def get_generation(self):
        """
        Return all the elements tagged as belonging to the current
        generation

        :return:
        """
        return [x for x in self.generation if self.current_generation in self.generation[x]]

    def get_generation_evaluated(self):
        """
        Return all the elements in the current generation that are also relaxed

        :return:
        """
        evaluated = self.population.evaluated
        return [x for x in self.population.actives if self.population.is_evaluated(x)]

    def pass_to_new_generation(self, imember):
        self.generation[imember].append(self.current_generation + 1)

    def save_generations(self):
        wf = open('generations.dat', 'w')
        for i in sorted(self.generation):
            wf.write('%15s %10d\n' % (i, len(self.generation[i])))
        wf.close()

    def replacing(self, imember, reason=None):
        self.population.disable(imember)
        new_member = self.population.add_random()
        # Replacing the element in the same generation
        ngeneration = self.generation[imember][-1]
        self.generation[new_member] = [ngeneration]
        self.generation[imember].remove(ngeneration)

    def print_status(self):
        print 'Generation', self.current_generation, sorted(self.get_generation())
        print 'Selection ', self.current_generation, sorted(self.get_generation_evaluated())
        print 'Actives   ', self.current_generation, sorted(self.population.actives)

    @abstractmethod
    def run_one_cycle(self):
        pass

    @abstractmethod
    def set_params(self, params):
        pass

    def run_all_cycles(self):
        """
        Execute the total number of cycles

        :return:
        """
        sleep_time = 30

        icycle = 0
        while True:
            print '\n GENERATION ', icycle

            log.debug('Enforcing the size of generation: %d' % self.generation_size)
            self.enforce_generation_size()

            while len(self.population.actives_evaluated) < self.fraction_evaluated * self.generation_size:
                log.debug("Population still not evaluated ")
                time.sleep(sleep_time)

            log.debug('Actives           : %10d' % len(self.population.actives))
            log.debug('Actives evaluated : %10d' % len(self.population.actives_evaluated))

            log.debug('Removing not evaluated...')
            for entry_id in self.population.actves_no_evaluated:
                self.replacing(entry_id, reason='no_evaluated')

            log.debug('Removing duplicates...')
            duplicates = self.population.check_duplicates()
            log.debug('%d' % len(duplicates))
            for entry_id in duplicates:
                self.replacing(entry_id, reason='duplicate')

            log.debug('Running one cycle...')
            self.run_one_cycle()

            best_member = self.population.ids_sorted(self.get_generation_evaluated())[0]
            print 'Best member is :', best_member
            print self.generation[best_member]
            print len(self.generation[best_member])
            if len(self.generation[best_member]) > self.stabilization_limit:
                break
            else:
                icycle += 1

        self.save_generations()

