from .searcher import Searcher
from pychemia import pcm_log


class ParticleSwarm(Searcher):
    def __init__(self, population, params=None, generation_size=32, stabilization_limit=10):
        """
        Implementation fo the Particle Swarm Optimization algorithm for global minimization

        :param population:
        :param params: (dict) Parameters to setup the Searcher. For PSO the parameters are:
            'beta0': A value in the range (0,1) that defines how close the mobile structure will be from to
             the target structure
            'alpha0': Factor of scale for the random movement
            'delta': How the random change decreases with time
        :param generation_size: (int)
        :param stabilization_limit: (int)
        :return:
        """
        # Initializing objects
        Searcher.__init__(self, population, generation_size, stabilization_limit)
        # Parameters
        self.beta0 = None
        self.alpha0 = None
        self.delta = None
        self.set_params(params)

    def set_params(self, params):
        # Default parameters
        self.beta0 = 0.9
        self.alpha0 = 0.3
        self.delta = 0.9
        if params is None:
            params = {}
        if 'beta0' in params:
            self.beta0 = params['beta0']
        if 'alpha0' in params:
            self.alpha0 = params['alpha0']
        if 'delta' in params:
            self.delta = params['delta']

    def get_params(self):
        return {'beta0': self.beta0,
                'alpha0': self.alpha0,
                'delta': self.delta}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.actives_in_generation)
        pcm_log.info('Size of selection : %d' % len(selection))

        # Minus sign because we are searching for minima
        intensity = self.population.get_values(selection)
        for entry_id in intensity:
            intensity[entry_id] *= -1

        moves = {}
        new_selection = {}
        for entry_id in selection:
            new_selection[entry_id] = None
            moves[entry_id] = 0

        # Move all the particles (Except the elite)
        # as the selection is sorted it means that the first one will no move
        pcm_log.debug('No Moving %d %s. Intensity: %7.3f' % (0, str(selection[0]), intensity[selection[0]]))
        for i in range(1, len(selection)):
            entry_id = selection[i]
            pcm_log.debug('Moving %d %s. Intensity: %7.3f' % (i, str(entry_id), intensity[entry_id]))

            new_selection[entry_id] = self.population.move(entry_id, selection[0], factor=self.beta0, in_place=False)
            factor = self.alpha0 * (self.delta ** self.current_generation)
            self.population.move_random(new_selection[entry_id], factor=factor, in_place=True)

        for entry_id in selection:
            if new_selection[entry_id] is not None:
                pcm_log.debug('[%s] Moved to: %s (%d moves)' %
                              (str(entry_id), new_selection[entry_id], moves[entry_id]))
                self.replace_by_other(entry_id, new_selection[entry_id],
                                      reason='Moved %d times' % moves[entry_id])
                self.population.enable(new_selection[entry_id])
            else:
                pcm_log.debug('[%s] Promoted to new generation' % str(entry_id))
                self.pass_to_new_generation(entry_id, reason='The best')
