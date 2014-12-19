__author__ = 'Guillermo Avendano-Franco'

import random
import time

from _genealogy import Genealogy


class HarmonySearch(Genealogy):

    def __init__(self, population, objective_function, evaluator, hmcr=0.9, par=0.9, top=2, tail=2,
                 stabilization_limit=5, fraction_evaluated=0.8, timeout_per_cycle=60):
        """
        Harmony Search Method:
        This searcher is the simplest one, does not require a metric space except for the evaluation
        of duplicates. It is more elaborated than a simple random search.
        Harmony Search only requires changers for existing members and generation of new random members.

        :param population: A subclass of the ABC Population
        :param objective_function: A subclass of the ABC ObjectiveFunction
        :param evaluator: A subclass of Evaluator
        :param hmcr: Harmony Memory Considering Rate, how often a value should be discarded in favor of new random
                    members, a value of 1 means never discard any member except those in tail (see below). A value
                    of 0 means always discard except those in tail (see below).
        :param par: Probability of being changed for those not discarded. A value of 1 means always changed, a value
                    of 0 means never changed.
        :param stabilization_limit: Number of generations survived by the best member
        :param fraction_evaluated: The ratio between the active members evaluated and the total number of members.
                    Once this ratio is fulfilled the generation is considered sufficient to start a cycle
        :param top: Number of members that are automatically promoted to the next generation (The best)
        :param tail: Number of members that are automatically discarded (The worst)
        :param timeout_per_cycle: Time in seconds before all the unevaluated members being discarded and new members
                    created
        """
        # Initializing variables
        self.population = population
        self.objective_function = objective_function
        self.evaluator = evaluator
        self.hmcr = None
        self.set_hmcr(hmcr)
        self.par = None
        self.set_par(par)
        self.stabilization_limit = stabilization_limit
        self.fraction_evaluated = fraction_evaluated
        self.top = top
        self.tail = tail
        self.timeout_per_cycle = timeout_per_cycle
        # Initializing objects
        self.objective_function.initialize(self.population)
        self.evaluator.initialize(self.population)
        Genealogy.__init__(self, self.population, self.evaluator)

    def set_hmcr(self, hmcr):
        """
        Check and set the harmony memory considering rate
        :param hmcr: A float value between 0 to 1
        """
        assert(1.0 >= hmcr >= 0.0)
        self.hmcr = hmcr

    def set_par(self, par):
        """
        Chack and set the probability of change
        :param par: A float value between 0 to 1
        """
        assert(1.0 >= par >= 0.0)
        self.par = par

    @property
    def harmony_memory_considering_rate(self):
        return self.hmcr

    def run_one_cycle(self):
        """
        Run one cycle for the Harmony Search Method
        """
        # Get a static selection of the values in the generation that are relaxed
        selection = self.get_generation_evaluated()
        sorted_selection = self.objective_function.ids_sorted(selection)
        print 'Size of selection : ', len(selection)
        print 'Size of actives   : ', len(self.population.actives)
        print 'Size of members   : ', len(self.population.members)
        print 'Size of generation   : ', len(self.generation)
        self.print_status()

        #Automatic promotion for active members that are not evaluated
        for imember in [x for x in self.population.actives if x not in sorted_selection]:
            print imember, " Active, not evaluated, promoted ",  self.population.member_str(imember)
            self.pass_to_new_generation(imember)

        #Automatic promotion for the top ranking members
        for imember in sorted_selection[:self.top]:
            print imember, " Good value, promoted ",  self.population.member_str(imember)
            self.pass_to_new_generation(imember)

        # Intermediate members, their fate depends on hmcr and par
        for imember in sorted_selection[self.top:-self.tail]:
            rnd = random.random()
            print imember, " In the middle (hmcr) %5.2f vs %5.2f)" % (rnd, self.hmcr)
            if rnd < self.hmcr:
                rnd = random.random()
                print imember, "  Promoted (par) %5.2f vs %5.2f)" % (rnd, self.par)
                if rnd < self.par:
                    self.population.disable(imember)
                    new_member = self.population.add_modified(imember)
                    self.generation[new_member] = [self.current_generation + 1]
                    print imember, '   Changed %s -> %s' % (imember, new_member)
                else:
                    print imember, '   Unchanged'
                    self.pass_to_new_generation(imember)
            else:
                print imember, '   Discarded '
                self.population.disable(imember)
                new_member = self.population.add_random()
                self.generation[new_member] = [self.current_generation + 1]

        for imember in sorted_selection[-self.tail:]:
            print imember, " Bad value, demoted ",  self.population.member_str(imember)
            self.population.disable(imember)
            new_member = self.population.add_random()
            self.generation[new_member] = [self.current_generation + 1]

        # Increase the current generation number
        self.current_generation += 1

    def run_all_cycles(self):
        """
        Execute the total number of cycles

        :return:
        """
        sleep_time = 2

        icycle = 0
        while True:
            print '\n GENERATION ', icycle

            if not self.evaluator.is_running:
                print 'Starting evaluator'
                self.evaluator.run()

            timeout = 0
            while True:
                if len(self.population.evaluated) > 0 and self.population.fraction_evaluated >= self.fraction_evaluated:
                    if self.population.fraction_evaluated < 1.0:
                        print "Some members lost: ", [x for x in self.population.actives
                                                      if x not in self.population.evaluated]
                        self.print_status()
                        for imember in self.population.active_no_evaluated():
                            print 'Removing: ', imember
                            self.replacing(imember)

                    print 'Checking duplicates...'
                    for imember in self.population.check_duplicates():
                        self.replacing(imember)
                    print 'Running one cycle...'
                    self.run_one_cycle()
                    break
                else:
                    print 'Fraction evaluated:', self.population.fraction_evaluated
                    timeout += sleep_time
                    time.sleep(sleep_time)
                    if timeout > self.timeout_per_cycle:
                        print 'Timeout for a single cycle, discarding unevaluated members'
                        for imember in self.population.active_no_evaluated():
                            print 'Removing: ', imember
                            self.print_status()
                            self.replacing(imember)
                    elif timeout > 2*self.timeout_per_cycle:
                        print 'Waiting too much, stopping now'
                        self.evaluator.stop()
                        return icycle
            best_member = self.objective_function.ids_sorted(self.get_generation_evaluated())[0]
            print 'Best member is :', best_member
            print self.generation[best_member]
            print len(self.generation[best_member])
            if len(self.generation[best_member]) > self.stabilization_limit:
                break
            else:
                icycle += 1

        self.save_generations()
        self.evaluator.stop()
