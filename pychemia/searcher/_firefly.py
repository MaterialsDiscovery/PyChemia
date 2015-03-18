__author__ = 'Guillermo Avendano-Franco'

import math
from pychemia import log
from _searcher import Searcher


class FireFly(Searcher):

    def __init__(self, population, params=None, fraction_evaluated=0.95, generation_size=32, stabilization_limit=10):
        """
        Implementation fo the Firefly algorithm for global minimization
        This searcher uses a metric to compute the attractiveness and the vector displacement
        to move one firefly in the direction of another one

        :param population:
        :param params: (dict) Parameters to setup the Searcher
        :param fraction_evaluated: (float)
        :param generation_size: (int)
        :param stabilization_limit: (int)
        :return:
        """
        # Mandatory objects
        self.population = population
        # Parameters
        self.gamma = None
        self.beta0 = None
        self.alpha0 = None
        self.delta = None
        self.multi_move = None
        self.set_params(params)
        # Constrains
        self.fraction_evaluated = fraction_evaluated
        self.generation_size = generation_size
        self.stabilization_limit = stabilization_limit
        # Initializing objects
        Searcher.__init__(self, self.population, fraction_evaluated, generation_size, stabilization_limit)

    def set_params(self, params):
        # Default parameters
        self.gamma = 0.1
        self.beta0 = 0.9
        self.alpha0 = 0.3
        self.delta = 0.9
        self.multi_move = False
        if 'gamma' in params:
            self.gamma = params['gamma']
        if 'beta0' in params:
            self.beta0 = params['beta0']
        if 'alpha0' in params:
            self.alpha0 = params['alpha0']
        if 'delta' in params:
            self.delta = params['delta']
        if 'multi_move' in params:
            self.multi_move = params['multi_move']

    def get_params(self):
        return {'gamma': self.gamma,
                'beta0': self.beta0,
                'alpha0': self.alpha0,
                'delta': self.delta,
                'multi_move': self.multi_move}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.population.actives_evaluated)

        # Minus sign because we are searching for minima
        intensity = self.population.get_values(selection)
        for entry_id in intensity:
            intensity[entry_id] *= -1

        moves = {}
        new_selection = {}
        for entry_id in selection:
            new_selection[entry_id] = None
            moves[entry_id] = 0

        # Move all the fireflies (Except the most brightness)
        # as the selection is sorted it means that the first one will no move
        log.debug('No Moving %d %s. Intensity: %7.3f' % (0, str(selection[0]), intensity[selection[0]]))
        for i in range(1, len(selection)):
            entry_id = selection[i]
            log.debug('Moving %d %s. Intensity: %7.3f' % (i, str(entry_id), intensity[entry_id]))

            # Moving in the direction of all the brighter fireflies
            if self.multi_move:

                for j in range(0, i):
                    entry_jd = selection[j]

                    distance = self.population.distance(entry_id, entry_jd)
                    beta = self.beta0*math.exp(-self.gamma * distance * distance)
                    # The variation of attractiveness \beta with the distance r
                    log.debug('[%s] Distance: %7.3f. Intensity: %7.3f. Atractiveness: %7.3f' %
                              (str(entry_jd), distance, intensity[entry_jd], beta))

                    if new_selection[entry_id] is None:
                        new_selection[entry_id] = self.population.move(entry_id, entry_jd, factor=beta, in_place=False)
                        self.population.move_random(new_selection[entry_id],
                                                    factor=self.alpha0*(self.delta**self.current_generation))
                    else:
                        self.population.move(new_selection[entry_id], entry_jd, in_place=True)
                        self.population.move_random(new_selection[entry_id],
                                                    factor=self.alpha0*(self.delta**self.current_generation))
                    moves[entry_id] += 1

            # Moving in the direction of the closets brighter firefly
            else:
                distances = [self.population.distance(entry_id, entry_jd) for entry_jd in selection[:i]]
                distance = min(distances)
                target = selection[distances.index(distance)]
                beta = self.beta0*math.exp(-self.gamma * distance * distance)
                # The variation of attractiveness \beta with the distance r
                log.debug('[%s] Distance: %7.3f. Intensity: %7.3f. Atractiveness: %7.3f' %
                          (str(entry_jd), distance, intensity[entry_jd], beta))

                new_selection[entry_id] = self.population.move(entry_id, target, factor=beta, in_place=False)
                self.population.move_random(new_selection[entry_id],
                                                factor=self.alpha0*(self.delta**self.current_generation))
                moves[entry_id] += 1

        for entry_id in selection:
            if new_selection[entry_id] is not None:
                log.debug('[%s] Moved to: %s ' % (str(entry_id), new_selection[entry_id]))
                self.replace_by_other(entry_id, new_selection[entry_id],
                                      reason='Moved %d times' % moves[entry_id])
                self.population.activate(new_selection[entry_id])
            else:
                log.debug('[%s] Promoted to new generation' % str(entry_id))
                self.pass_to_new_generation(entry_id, reason='The best')
