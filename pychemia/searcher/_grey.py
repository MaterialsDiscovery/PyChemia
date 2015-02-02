__author__ = 'Guillermo Avendano-Franco'

import time

from _genealogy import Genealogy


class GreyWolf(Genealogy):
    def __init__(self, population, evaluator, objective_function, stabilization_limit=10, fraction_evaluated=0.8,
                 param_a=1, param_c=1):
        self.population = population
        self.evaluator = evaluator
        self.objective_function = objective_function
        self.stabilization_limit = stabilization_limit
        self.fraction_evaluated = fraction_evaluated
        Genealogy.__init__(self, self.population, self.evaluator)
        self.param_a = None
        self.set_param_a(param_a)
        self.param_c = None
        self.set_param_c(param_c)

    def set_param_a(self, param_a):
        self.param_a = param_a

    def set_param_c(self, param_c):
        self.param_c = param_c

    def run_one_cycle(self):
        pass

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
                if len(self.population.evaluated) > 0 and self.population.fraction_evaluated > self.fraction_evaluated:
                    if self.population.fraction_evaluated < 1.0:
                        print "Some members lost"
                        self.print_status()
                        for imember in self.population.active_no_evaluated:
                            print 'Removing: ', imember
                            self.replacing(imember)

                    for imember in self.population.check_duplicates():
                        self.replacing(imember)
                    self.run_one_cycle()
                    break
                else:
                    timeout += sleep_time
                    time.sleep(sleep_time)
                    if timeout > self.timeout_per_cycle:
                        print 'Timeout for a single cycle, discarding unevaluated members'
                        for imember in self.population.active_no_evaluated:
                            print 'Removing: ', imember
                            self.print_status()
                            self.replacing(imember)
                    elif timeout > 2 * self.timeout_per_cycle:
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
