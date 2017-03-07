import math
import numpy as np

from .searcher import Searcher
from pychemia import pcm_log


class FireFly(Searcher):
    def __init__(self, population, params=None, generation_size=32,
                 stabilization_limit=10, target_value=None):
        """
        Implementation fo the Firefly algorithm for global minimization
        This searcher uses a metric to compute the attractiveness and the vector displacement
        to move one firefly in the direction of another one

        :param population:
        :param params: (dict) Parameters to setup the Searcher. For Firefly the parameters are:
            'beta0': A value in the range (0,1) that defines how close the mobile structure will be from to
             the target structure
            'gamma': Defines the exponent of decreasing for the movement of the mobile structure
            'alpha0': Factor of scale for the random movement
            'delta': How the random change decreases with time
            'multi_move': Boolean to express if the fireflies moves following all the  other brighter ones or just
            the closest brigter firefly
        :param generation_size: (int)
        :param stabilization_limit: (int)
        :return:
        """
        # Initializing objects
        Searcher.__init__(self, population, generation_size, stabilization_limit, target_value)
        # Parameters
        self.gamma = None
        self.beta0 = None
        self.alpha0 = None
        self.delta = None
        self.multi_move = None
        self.set_params(params)

    def set_params(self, params):
        # Default parameters
        self.gamma = 0.1
        self.beta0 = 0.9
        self.alpha0 = 0.3
        self.delta = 0.9
        self.multi_move = True
        if params is None:
            params = {}
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
        selection = self.population.ids_sorted(self.actives_in_generation)
        pcm_log.info(' Size of selection : %d' % len(selection))

        # For statistical purposes
        distances = []
        intensities = []
        atractiveness = []

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
        pcm_log.debug('No Moving %d %s. Intensity: %7.3f' % (0, str(selection[0]), intensity[selection[0]]))
        for i in range(1, len(selection)):
            entry_id = selection[i]
            pcm_log.debug('Moving %d %s. Intensity: %7.3f' % (i, str(entry_id), intensity[entry_id]))

            # Moving in the direction of all the brighter fireflies
            if self.multi_move:

                for j in list(range(0, i))[::-1]:
                    entry_jd = selection[j]

                    distance = self.population.distance(entry_id, entry_jd)
                    beta = self.beta0 * math.exp(-self.gamma * distance * distance)
#                    The variation of attractiveness \beta with the distance r
                    pcm_log.debug('[%s] Distance: %7.3f. Intensity: %7.3f. Atractiveness: %7.3f' % (str(entry_jd),
                                                                                                    distance,
                                                                                                    intensity[entry_jd],
                                                                                                    beta))
                    # Collecting Statistics
                    distances.append(distance)
                    intensities.append(intensity[entry_jd])
                    atractiveness.append(beta)

                    if new_selection[entry_id] is None:
                        new_selection[entry_id] = self.population.move(entry_id, entry_jd, factor=beta, in_place=False)
                        if self.alpha0 > 0:
                            factor = self.alpha0 * (self.delta ** self.current_generation)
                            self.population.move_random(new_selection[entry_id], factor=factor, in_place=True)
                    else:
                        self.population.move(new_selection[entry_id], entry_jd, factor=beta, in_place=True)
                        if self.alpha0 > 0:
                            factor = self.alpha0 * (self.delta ** self.current_generation)
                            self.population.move_random(new_selection[entry_id], factor=factor, in_place=True)
#                    print(new_selection)
                    moves[entry_id] += 1

            # Moving in the direction of the closets brighter firefly
            else:
                distances = [self.population.distance(entry_id, entry_jd) for entry_jd in selection[:i]]
                distance = min(distances)
                target = selection[distances.index(distance)]
                beta = self.beta0 * math.exp(-self.gamma * distance * distance)
                # The variation of attractiveness \beta with the distance r
                pcm_log.debug('[%s] Distance: %7.3f. Intensity: %7.3f. Atractiveness: %7.3f' %
                              (str(entry_jd), distance, intensity[entry_jd], beta))

                new_selection[entry_id] = self.population.move(entry_id, target, factor=beta, in_place=False)
                factor = self.alpha0 * (self.delta ** self.current_generation)
                self.population.move_random(new_selection[entry_id], factor=factor, in_place=True)
                moves[entry_id] += 1

        if len(distances) > 0:
            pcm_log.info('+----------------+--------------+-------------+-------------+')
            pcm_log.info('+                |    Minimum   |   Maximum   |   Average   |')
            pcm_log.info('+----------------+--------------+-------------+-------------+')
            pcm_log.info('+ Distances      |    %7.2f   |   %7.2f   |   %7.2f   |' % (np.min(distances),
                                                                                      np.max(distances),
                                                                                      np.average(distances)))
            pcm_log.info('+ Intensities    |    %7.2f   |   %7.2f   |   %7.2f   |' % (np.min(intensities),
                                                                                      np.max(intensities),
                                                                                      np.average(intensities)))
            pcm_log.info('+ Attractiveness |    %7.2f   |   %7.2f   |   %7.2f   |' % (np.min(atractiveness),
                                                                                      np.max(atractiveness),
                                                                                      np.average(atractiveness)))
            pcm_log.info('+----------------+--------------+-------------+-------------+')

        for entry_id in selection:
            if new_selection[entry_id] is not None:
                self.replace_by_other(entry_id, new_selection[entry_id],
                                      reason='Moved %d times' % moves[entry_id])
                self.population.enable(new_selection[entry_id])
            else:
                self.pass_to_new_generation(entry_id, reason='The best')
