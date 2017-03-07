from numpy import sum
from .searcher import Searcher
from pychemia import pcm_log

"""
Implementation of Genetic Algorithms
This is the abstract layer
No references to Structural Search should be appear
This module should be general enough even for stock market prediction!
"""

__author__ = 'Guillermo Avendano-Franco'


class GeneticAlgorithm(Searcher):
    def __init__(self, population, params=None, generation_size=32, stabilization_limit=10):
        Searcher.__init__(self, population, generation_size, stabilization_limit)
        self.nelite = None
        self.mutation_prob = None
        self.cross_prob = None
        self.crossing_sets = None
        self.set_params(params)

    def set_params(self, params):
        self.nelite = 2
        self.mutation_prob = 0.5
        self.cross_prob = 0.5
        if self.generation_size == 32:
            self.crossing_sets = [10, 5]
        elif self.generation_size == 16:
            self.crossing_sets = [5, 2]
        elif self.generation_size == 8:
            self.crossing_sets = [2, 1]

        if params is None:
            params = {}
        if 'nelite' in params:
            self.nelite = params['nelite']
        if 'mutation_prob' in params:
            self.mutation_prob = params['mutation_prob']
        if 'cross_prob' in params:
            self.cross_prob = params['cross_prob']
        if 'crossing_sets' in params:
            self.crossing_sets = params['crossing_sets']

        self.nelite = int(self.generation_size - 2 * sum(self.crossing_sets))
        if self.nelite <= 0:
            raise ValueError('Enter explicit crossing_sets for your generation size')

    def get_params(self):
        return {'nelite': self.nelite, 'mutation_prob': self.mutation_prob, 'cross_prob': self.cross_prob,
                'crossing_sets': self.crossing_sets}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.actives_in_generation)
        pcm_log.info(' Size of selection : %d' % len(selection))

        # Minus sign because we are searching for minima
        intensity = self.population.get_values(selection)
        for entry_id in intensity:
            intensity[entry_id] *= -1

        discarded_index = self.nelite + int(sum(self.crossing_sets))

        for i in range(min(self.nelite, len(selection))):
            entry_id = selection[i]
            pcm_log.debug('[%s] Promoted to new generation' % str(entry_id))
            self.pass_to_new_generation(entry_id, reason='Elite')

        jump = 0
        for i in range(min(len(self.crossing_sets), len(selection))):
            for j in range(self.crossing_sets[i]):
                entry_id = selection[i]

                if self.nelite + j + jump < len(selection):
                    entry_jd = selection[self.nelite + j + jump]
                    new_entry_id, new_entry_jd = self.population.cross([entry_id, entry_jd])
                    msg = 'Replace candidates %d and %d by crossing %d with %d'
                    pcm_log.info(msg % (self.nelite + j + jump, discarded_index, self.nelite + j + jump, i))

                    pcm_log.debug('[%s] Moved to: %s' % (entry_id, new_entry_id))
                    self.replace_by_other(entry_jd, new_entry_id,
                                          reason='Cross between %s and %s' % (entry_id, entry_jd))
                    self.population.enable(new_entry_id)
                    if discarded_index < len(selection):
                        entry_kd = selection[discarded_index]
                        pcm_log.debug('[%s] Moved to: %s' % (entry_kd, new_entry_jd))
                        self.replace_by_other(entry_kd, new_entry_jd,
                                              reason='Cross between %s and %s' % (entry_id, entry_jd))
                        self.population.enable(new_entry_jd)
                    else:
                        pcm_log.debug('Candidate %s will not be activated' % new_entry_jd)

                discarded_index += 1
            jump += self.crossing_sets[i]
