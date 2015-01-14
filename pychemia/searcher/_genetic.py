"""
Implementation of Genetic Algorithms
This is the abstract layer
No references to Structural Search should be appear
This module should be general enough even for stock market prediction!
"""

from _genealogy import Genealogy
from _searcher import Searcher

class GeneticAlgorithm(Genealogy, Searcher):

    def __init__(self, population, evaluator, objective_function, params, stabilization_limit=10,
                 fraction_evaluated=0.8):
        self.population = population
        self.evaluator = evaluator
        Genealogy.__init__(self, self.population, self.evaluator)

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
                        for imember in self.population.active_no_evaluated():
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
