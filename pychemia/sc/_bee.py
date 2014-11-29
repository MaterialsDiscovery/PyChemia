__author__ = 'Guillermo Avendano-Franco'

from _metaheuristics import MetaHeuristics


class BeeAlgorithm(MetaHeuristics):

    def __init__(self, population, objective_function, relax_population, ne, nre, nb, nrb, ns, number_generations,
                 ratio_relax=0.8):
        """
        Implementation fo the Firefly algorithm for global minimization

        :param population:
        :param objective_function:
        :param relax_population:
        :param gamma:
        :param number_generations:
        :param ratio_relax:
        :return:
        """
        self.population = population
        self.objective_function = objective_function
        self.relax_population = relax_population
        self.ne = None   # Number of elite scout bees
        self.nre = None  # Number of elite foragers
        self.nb = None   # Number of best scout bees (excluding elites)
        self.nrb = None  # Number of best foragers
        self.ns = None   # Number of scouts (include elites,  best and others)
        self.n = None    # Number of bees (size of colony n=ne*nre + nb*nrb + ns )
        self.set_parameters(ne, nre, nb, nrb, ns)
        self.number_generations = number_generations
        self.ratio_relax = ratio_relax
        self.objective_function.initialize(self.population)
        self.relax_population.initialize(self.population)
        MetaHeuristics.__init__(self, self.population, self.relax_population)
        self.scouts_elite = None
        self.scouts_best = None
        self.scouts_others = None
        self.foragers = {}

    def set_parameters(self, ne, nre, nb, nrb, ns):
        """
        Set all the parameters required by the bee algorithm,
        Check right values for all the parameters

        :param ne:
        :param nre:
        :param nb:
        :param nrb:
        :param ns:
        """
        assert(nre > 0)
        assert(nrb > 0)
        assert(ne > 0)
        assert(nb > 0)
        assert(nb > ne)
        assert(ns > nb)
        self.nb = nb
        self.nrb = nrb
        self.ne = ne
        self.nre = nre
        self.ns = ns
        self.n = ne*nre+(nb-ne)*nrb+ns

    def run_one_cycle(self):

        # Get a static selection of the values in the generation that are relaxed
        selection = self.get_generation_relaxed()

        if self.current_generation == 0:
            # During the first iteration all the selection are scouts
            # Lets choose from there the elite, the best and the actual
            # scouts for the next iterations
            if len(selection) >= self.ns:
                members = self.objective_function.sort_by_fitness(selection)
                self.scouts_elite = members[:self.ne]
                self.scouts_best = members[self.ne:self.ne+self.nb]
                self.scouts_others = members[self.ne+self.nb:self.ns-self.ne+self.nb]
                # Disable extra scouts (from an initial population)
                for imember in members[self.ns:]:
                    self.population.disable(imember)

                # Add nre foragers around each elite scout
                for imember in self.scouts_elite:
                    for i in range(self.nre):
                        forager = self.population.add_modified(imember)
                        if imember not in self.foragers:
                            self.foragers['imember'] = [forager]
                        else:
                            self.foragers['imember'].append(forager)

                # Add nrb foragers around each best scout
                for imember in self.scouts_best:
                    for i in range(self.nrb):
                        forager = self.population.add_modified(imember)
                        if imember not in self.foragers:
                            self.foragers['imember'] = [forager]
                        else:
                            self.foragers['imember'].append(forager)
            else:
                # Number of scouts insufficient, increase their number with more
                # random members
                for i in range(selection, self.ns):
                    self.population.add_random()
        else:
            # New iterations look into each patch and see witch bees are on each patch
            for i in self.scouts_elite:
                new_elite = i
                dead_bees = 0
                for j in self.foragers[i].copy():
                    if j not in selection:
                        # Consider a dead bee
                        self.foragers[i].remove(j)
                        dead_bees += 1
                members = self.objective_function.sort_by_fitness(self.foragers[i])