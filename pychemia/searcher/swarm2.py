import math

from .searcher import Searcher
from pychemia import pcm_log


class ParticleSwarm(Searcher):
    def __init__(self, population, params=None, generation_size=32, stabilization_limit=10):
        """
        Implementation fo the Firefly algorithm for global minimization
        This searcher uses a metric to compute the attractiveness and the vector displacement
        to move one firefly in the direction of another one

        :param population:
        :param params: (dict) Parameters to setup the Searcher
        :param generation_size: (int)
        :param stabilization_limit: (int)
        :return:
        """
        # Mandatory objects
        self.population = population
        # Parameters
        self.gamma = None
        self.elites = None
        self.set_params(params)
        # Constrains
        self.generation_size = generation_size
        self.stabilization_limit = stabilization_limit
        # Initializing objects
        Searcher.__init__(self, self.population, generation_size, stabilization_limit)

    def set_params(self, params):
        if params is None:
            self.gamma = 0.1
            self.elites = 3
        else:
            assert ('gamma' in params)
            assert (params['gamma'] >= 0.0)
            self.gamma = params['gamma']
            if 'elites' in params:
                self.elites = params['elites']

    def get_params(self):
        return {'gamma': self.gamma, 'elites': self.elites}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.actives_in_generation)

        # Minus sign because we are searching for minima
        intensity = self.population.get_values(selection)
        for entry_id in intensity:
            intensity[entry_id] *= -1

        moves = {}
        new_selection = {}
        for entry_id in selection:
            new_selection[entry_id] = None

        # Move all the fireflies (Except the most brightness)
        # as the selection is sorted it means that the first one will no move
        pcm_log.debug('No Moving %d %s. Intensity: %7.3f' % (0, str(selection[0]), intensity[selection[0]]))

        # The best
        elites = selection[:self.elites]

        for i in range(self.elites, len(selection)):
            entry_id = selection[i]
            pcm_log.debug('Moving %d %s. Intensity: %7.3f' % (i, str(entry_id), intensity[entry_id]))

            distances = [self.population.distance(entry_id, entry_jd) for entry_jd in elites]
            target = elites[distances.index(min(distances))]
            distance = min(distances)
            atractiveness = math.exp(-self.gamma * distance) * intensity[target]

            pcm_log.debug('[%s] Distance: %7.3f. Intensity: %7.3f. Atractiveness: %7.3f' % (str(target),
                                                                                            distance,
                                                                                            intensity[target],
                                                                                            atractiveness))

            if intensity[entry_id] < atractiveness:
                new_selection[entry_id] = self.population.move(entry_id, target, in_place=False)

        for entry_id in selection:
            pcm_log.debug('Deciding fate for firefly: %s' % str(entry_id))
            if new_selection[entry_id] is not None:
                pcm_log.debug('Moved to a new location %s ' % str(entry_id))
                self.replace_by_other(entry_id, new_selection[entry_id], reason=None)
            else:
                pcm_log.debug('Promoted to new generation ')
                self.pass_to_new_generation(entry_id, reason='No other firefly is more attractive')
