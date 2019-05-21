
from .searcher import Searcher
from pychemia import pcm_log


class BeeAlgorithm(Searcher):
    def __init__(self, population, params=None, generation_size=32, stabilization_limit=10):
        """
        Implementation fo the Firefly algorithm for global minimization

        :param population:
        :param stabilization_limit:
        :return:
        """
        Searcher.__init__(self, population, generation_size, stabilization_limit)
        # Parameters
        self.ne = None  # Number of elite scout bees
        self.nre = None  # Number of elite foragers
        self.nb = None  # Number of best scout bees (excluding elites)
        self.nrb = None  # Number of best foragers
        self.ns = None  # Number of scouts (include elites,  best and others)
        self.n = None  # Number of bees (size of colony n=ne*nre + nb*nrb + ns )
        self.set_params(params)
        self.scouts_elite = None
        self.scouts_best = None
        self.scouts_others = None
        self.foragers = {}
        self.delta_change = 1

    def set_params(self, params):
        """
        Set all the parameters required by the bee algorithm,
        Check right values for all the parameters

        :param params:
        """
        self.n = self.generation_size  # Number of bees (size of colony n=ne*nre + nb*nrb + ns )
        self.ne = 3  # Number of elite scout bees
        self.nre = 3  # Number of elite foragers
        self.nb = 2  # Number of best scout bees (excluding elites)
        self.nrb = 2  # Number of best foragers
        # Number of scouts (include elites,  best and others)
        self.ns = self.n - self.ne * self.nre - self.nb * self.nrb
        if params is None:
            params = {}
        if 'nb' in params:
            self.nb = params['nb']
        if 'nrb' in params:
            self.nrb = params['nrb']
        if 'ne' in params:
            self.ne = params['ne']
        if 'nre' in params:
            self.nre = params['nre']
        assert (self.ns < self.n)

    def get_params(self):
        return {'nb': self.nb, 'ne': self.ne, 'nrb': self.nrb, 'nre': self.nre}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.actives_in_generation)
        pcm_log.info('Size of selection : %d' % len(selection))

        if self.scouts_elite is None:
            print('First run', len(selection), self.ns)
            # During the first iteration all the selection are scouts
            # Lets choose from there the elite, the best and the actual
            # scouts for the next iterations
            if len(selection) >= self.ns:
                self.scouts_elite = list(selection[:self.ne])
                self.scouts_best = list(selection[self.ne:self.ne + self.nb])
                self.scouts_others = list(selection[self.ne + self.nb:self.ns - self.ne + self.nb])
                # Disable extra scouts (from an initial population)
                for entry_id in selection[self.ns:]:
                    self.population.disable(entry_id)
                self.print_status()

                # Add nre foragers around each elite scout
                for entry_id in self.scouts_elite:
                    self.pass_to_new_generation(entry_id, reason='Elite')
                    self.create_foragers(entry_id, self.nre)
                self.print_status()

                # Add nrb foragers around each best scout
                for entry_id in self.scouts_best:
                    self.pass_to_new_generation(entry_id, reason='Best')
                    self.create_foragers(entry_id, self.nrb)
                self.print_status()

                for u in range(self.generation_size - len(self.get_generation(self.current_generation + 1))):
                    ident, origin = self.population.add_random()
                    self.population.disable(ident)
                    self.generation[ident] = [self.current_generation + 1]
                    self.scouts_others.append(ident)
                    print('Added to other scouts', ident)
                self.print_status()

            else:
                pass
                # Number of scouts insufficient, increase their number with more
                # random members
                # for i in range(selection, self.ns):
                #    entry_id=self.population.add_random()
                #    print 'Population raised:',entry_id
        else:
            # New iterations look into each patch and see witch bees are on each patch
            print('Foragers %d' % len(self.foragers), self.foragers)
            print('Elite    %d' % len(self.scouts_elite), self.scouts_elite)
            print('Best     %d' % len(self.scouts_best), self.scouts_best)
            print('Others   %d' % len(self.scouts_others), self.scouts_others)

            self.process_scouts(self.scouts_elite, selection)
            self.process_scouts(self.scouts_best, selection)

            dead_bees = 0
            for j in list(self.scouts_others):
                if j not in selection:
                    # Consider a dead bee
                    self.scouts_others.remove(j)
                    dead_bees += 1
            scouts_others_alive = self.population.ids_sorted(self.scouts_others)
            if len(scouts_others_alive) > 0:
                best_scout = scouts_others_alive[0]
                worst_best = self.population.ids_sorted(self.scouts_best)[-1]
                if self.population.value(best_scout) < self.population.value(worst_best):
                    self.scouts_best.remove(worst_best)
                    self.scouts_best.append(best_scout)
                    self.population.disable(worst_best)

            # Add nre foragers around each elite scout
            for entry_id in self.scouts_elite:
                self.pass_to_new_generation(entry_id, reason='Elite')
                self.foragers[entry_id] = []
                self.create_foragers(entry_id, self.nre)
            self.print_status()

            # Add nrb foragers around each best scout
            for entry_id in self.scouts_best:
                self.pass_to_new_generation(entry_id, reason='Best')
                self.foragers[entry_id] = []
                self.create_foragers(entry_id, self.nrb)
            self.print_status()

            for i in self.scouts_others:
                self.population.disable(i)
            self.print_status()

            # for i in range(self.ns - self.ne - self.nb):
            for u in range(self.generation_size - len(self.get_generation(self.current_generation + 1))):
                ident, origin = self.population.add_random()
                self.population.disable(ident)
                self.generation[ident] = [self.current_generation + 1]
                self.scouts_others.append(ident)
                print('Added to other scouts', ident)

            print('After run_one:')
            self.print_status()

    def create_foragers(self, scout, n):
        for i in range(n):
            forager = self.population.move_random(scout, factor=self.delta_change, in_place=False, kind='move')
            print('For scout:', scout, ' new forager: ', forager)
            self.population.disable(forager)
            self.generation[forager] = [self.current_generation + 1]
            if scout not in self.foragers:
                self.foragers[scout] = [forager]
            else:
                self.foragers[scout].append(forager)
                # change={'to': scout}
                # self.write_change(forager, change)

    def process_scouts(self, scouts, selection):
        # Removing dead bees from log
        for ibee in list(scouts):
            print(ibee, len(self.foragers[ibee]))
            dead_bees = 0
            for j in list(self.foragers[ibee]):
                if j not in selection:
                    # Consider a dead bee
                    self.foragers[ibee].remove(j)
                    dead_bees += 1

        for ibee in list(scouts):
            foragers_alive = self.population.ids_sorted(self.foragers[ibee])
            print('Foragers that survived:', foragers_alive)
            if len(foragers_alive) > 0:
                best_forager = foragers_alive[0]
                if self.population.value(best_forager) < self.population.value(ibee):
                    # We found a better elite
                    scouts.remove(ibee)
                    scouts.append(best_forager)
                    self.population.disable(ibee)
                    for i in foragers_alive[1:]:
                        self.population.disable(i)

                else:
                    for i in foragers_alive[:]:
                        self.population.disable(i)
