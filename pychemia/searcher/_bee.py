__author__ = 'Guillermo Avendano-Franco'

from _searcher import Searcher


class BeeAlgorithm(Searcher):

    def __init__(self, population, params, fraction_evaluated=0.8, generation_size=32, stabilization_limit=10):
        """
        Implementation fo the Firefly algorithm for global minimization

        :param population:
        :param stabilization_limit:
        :param fraction_evaluated:
        :return:
        """
        self.population = population
        self.ne = None   # Number of elite scout bees
        self.nre = None  # Number of elite foragers
        self.nb = None   # Number of best scout bees (excluding elites)
        self.nrb = None  # Number of best foragers
        self.ns = None   # Number of scouts (include elites,  best and others)
        self.n = None    # Number of bees (size of colony n=ne*nre + nb*nrb + ns )
        self.set_params(params)
        Searcher.__init__(self, self.population, fraction_evaluated, generation_size, stabilization_limit)
        self.scouts_elite = None
        self.scouts_best = None
        self.scouts_others = None
        self.foragers = {}

    def set_params(self, params):
        """
        Set all the parameters required by the bee algorithm,
        Check right values for all the parameters

        :param params:
        """
        assert(params['nre'] > 0)
        assert(params['nrb'] > 0)
        assert(params['ne'] > 0)
        assert(params['nb'] > 0)
        assert(params['nb'] > params['ne'])
        assert(params['ns'] > params['nb'])
        self.nb = params['nb']
        self.nrb = params['nrb']
        self.ne = params['ne']
        self.nre = params['nre']
        self.ns = params['ns']
        self.n = self.ne*self.nre+(self.nb-self.ne)*self.nrb+self.ns

    def run_one_cycle(self):

        # Get a static selection of the values in the generation that are relaxed
        selection = self.get_generation_evaluated()

        if self.scouts_elite is None:
            # During the first iteration all the selection are scouts
            # Lets choose from there the elite, the best and the actual
            # scouts for the next iterations
            if len(selection) >= self.ns:
                members = self.population.ids_sorted(selection)
                self.scouts_elite = members[:self.ne]
                self.scouts_best = members[self.ne:self.ne+self.nb]
                self.scouts_others = members[self.ne+self.nb:self.ns-self.ne+self.nb]
                # Disable extra scouts (from an initial population)
                for imember in members[self.ns:]:
                    self.population.disable(imember)

                # Add nre foragers around each elite scout
                for imember in self.scouts_elite:
                    self.create_foragers(imember, self.nre)

                # Add nrb foragers around each best scout
                for imember in self.scouts_best:
                    self.create_foragers(imember, self.nrb)
            else:
                # Number of scouts insufficient, increase their number with more
                # random members
                for i in range(selection, self.ns):
                    self.population.add_random()
        else:
            # New iterations look into each patch and see witch bees are on each patch
            for ielite in self.scouts_elite.copy():
                dead_bees = 0
                for j in self.foragers[ielite].copy():
                    if j not in selection:
                        # Consider a dead bee
                        self.foragers[ielite].remove(j)
                        dead_bees += 1
                foragers_alive = self.population.ids_sorted(self.foragers[ielite])
                best_forager = foragers_alive[0]
                if self.population.value(best_forager) < self.population.value(ielite):
                    # We found a better elite
                    self.scouts_elite.remove(ielite)
                    self.scouts_elite.append(best_forager)
                    self.population.disable(ielite)
                    for i in foragers_alive[1:]:
                        self.population.disable(i)
                else:
                    for i in foragers_alive[:]:
                        self.population.disable(i)
            for ibest in self.scouts_best.copy():
                dead_bees = 0
                for j in self.foragers[ibest].copy():
                    if j not in selection:
                        # Consider a dead bee
                        self.foragers[ibest].remove(j)
                        dead_bees += 1
                foragers_alive = self.population.ids_sorted(self.foragers[ibest])
                best_forager = foragers_alive[0]
                if self.population.value(best_forager) < self.population.value(ibest):
                    # We found a better elite
                    self.scouts_best.remove(ibest)
                    self.scouts_best.append(best_forager)
                    self.population.disable(ibest)
                    for i in foragers_alive[1:]:
                        self.population.disable(i)
                else:
                    for i in foragers_alive[:]:
                        self.population.disable(i)

            dead_bees = 0
            for j in self.scouts_others.copy():
                if j not in selection:
                    # Consider a dead bee
                    self.scouts_others.remove(j)
                    dead_bees += 1
            scouts_others_alive = self.population.ids_sorted(self.scouts_others)
            best_scout = scouts_others_alive[0]
            worst_best = self.population.ids_sorted(self.scouts_best)[-1]
            if self.population.value(best_scout) < self.population.value(worst_best):
                self.scouts_best.remove(worst_best)
                self.scouts_best.append(best_scout)
                self.population.disable(worst_best)

            # Add nre foragers around each elite scout
            for imember in self.scouts_elite:
                self.create_foragers(imember, self.nre)

            # Add nrb foragers around each best scout
            for imember in self.scouts_best:
                self.create_foragers(imember, self.nrb)

            for i in self.scouts_others:
                self.population.disable(i)

            for i in range(self.ns - self.ne - self.nb):
                ident = self.population.add_random()
                self.scouts_others.append(ident)

    def create_foragers(self, scout, n):
        for i in range(n):
            forager = self.population.add_modified(scout)
            if scout not in self.foragers:
                self.foragers[scout] = [forager]
            else:
                self.foragers[scout].append(forager)
