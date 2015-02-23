__author__ = 'Guillermo Avendano-Franco'

import random
from pychemia import log
from _searcher import Searcher


class HarmonySearch(Searcher):
    def __init__(self, population, params, fraction_evaluated=0.8, generation_size=32, stabilization_limit=10):
        """
        Harmony Search Method:
        This searcher is the simplest one, does not require a metric space except for the evaluation
        of duplicates. It is more elaborated than a simple random search.
        Harmony Search only requires changers for existing members and generation of new random members.

        :param population: A subclass of the ABC Population
        :param params: (dict) Parameters to setup the Searcher
                    hmcr: Harmony Memory Considering Rate, how often a value should be discarded in favor of new random
                    members, a value of 1 means never discard any member except those in tail (see below). A value
                    of 0 means always discard except those in tail (see below).
                    par: Probability of being changed for those not discarded. A value of 1 means always changed, a value
                    of 0 means never changed.
                    top: Number of members that are automatically promoted to the next generation (The best)
                    tail: Number of members that are automatically discarded (The worst)
        :param stabilization_limit: Number of generations survived by the best member
        :param fraction_evaluated: The ratio between the active members evaluated and the total number of members.
                    Once this ratio is fulfilled the generation is considered sufficient to start a cycle
        """
        # Mandatory objects
        self.population = population
        # Parameters
        self.hmcr = None  # harmony_memory_considering_rate
        self.par = None
        self.top = None
        self.tail = None
        self.set_params(params)
        # Initializing objects 
        Searcher.__init__(self, self.population, fraction_evaluated, generation_size, stabilization_limit)

    def set_params(self, params):
        if params is None:
            params = {'hmcr': 0.9, 'par': 0.9, 'top': 2, 'tail': 2}
        elif 'hmcr' not in params:
            params['hmcr'] = 0.9
        elif 'par' not in params:
            params['par'] = 0.9
        elif 'top' not in params:
            params['top'] = 2
        elif 'tail' not in params:
            params['tail'] = 2

        assert (1.0 >= params['hmcr'] >= 0.0)
        self.hmcr = params['hmcr']
        assert (1.0 >= params['par'] >= 0.0)
        self.par = params['par']

        self.top = params['top']
        self.tail = params['tail']

    def run_one_cycle(self):
        """
        Run one cycle for the Harmony Search Method
        """
        # Get a static selection of the values in the generation that are evaluated
        selection = self.get_generation_evaluated()
        sorted_selection = self.population.ids_sorted(selection)
        print 'Size of selection : ', len(selection)
        print 'Size of actives   : ', len(self.population.actives)
        print 'Size of members   : ', len(self.population.members)
        print 'Size of generation   : ', len(self.generation)
        self.print_status()

        # Automatic promotion for active members that are not evaluated
        for entry_id in [x for x in self.population.actives if x not in sorted_selection]:
            print entry_id, " Active, not evaluated, promoted "
            self.pass_to_new_generation(entry_id)

        # Automatic promotion for the top ranking members
        for entry_id in sorted_selection[:self.top]:
            print entry_id, " Good value, promoted "
            self.pass_to_new_generation(entry_id)

        # Intermediate members, their fate depends on hmcr and par
        for entry_id in sorted_selection[self.top:-self.tail]:
            rnd = random.random()
            print entry_id, " In the middle (hmcr) %5.2f vs %5.2f)" % (rnd, self.hmcr)
            if rnd < self.hmcr:
                rnd = random.random()
                print entry_id, "  Promoted (par) %5.2f vs %5.2f)" % (rnd, self.par)
                if rnd < self.par:
                    self.population.disable(entry_id)
                    new_member = self.population.add_modified(entry_id)
                    self.generation[new_member] = [self.current_generation + 1]
                    print entry_id, '   Changed %s -> %s' % (entry_id, new_member)
                else:
                    print entry_id, '   Unchanged'
                    self.pass_to_new_generation(entry_id)
            else:
                print entry_id, '   Discarded '
                self.population.disable(entry_id)
                new_member = self.population.add_random()
                self.generation[new_member] = [self.current_generation + 1]

        for entry_id in sorted_selection[-self.tail:]:
            print entry_id, " Bad value, demoted "
            self.population.disable(entry_id)
            new_member = self.population.add_random()
            self.generation[new_member] = [self.current_generation + 1]

        # Increase the current generation number
        self.current_generation += 1
