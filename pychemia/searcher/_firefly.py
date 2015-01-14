__author__ = 'Guillermo Avendano-Franco'

import math
import time

from _genealogy import Genealogy
from _searcher import Searcher

class FireFly(Genealogy, Searcher):

    def __init__(self, population, evaluator, objective_function, params, stabilization_limit=10,
                 fraction_evaluated=0.8):
        """
        Implementation fo the Firefly algorithm for global minimization
        This searcher uses a metric to compute the attractiveness and the vector displacement
        to move one firefly in the direction of another one

        :param population:
        :param objective_function:
        :param evaluator:
        :param gamma:
        :param number_generations:
        :param fraction_evaluated:
        :return:
        """
        self.population = population
        self.objective_function = objective_function
        self.evaluator = evaluator
        self.gamma = None
        self.set_params(params)
        self.stabilization_limit = stabilization_limit
        self.fraction_evaluated = fraction_evaluated

        self.objective_function.initialize(self.population)
        self.evaluator.initialize(self.population)
        Genealogy.__init__(self, self.population, self.evaluator)

    def set_params(self, params):
        assert('gamma' in params)
        assert(params['gamma'] >= 0.0)
        self.gamma = params['gamma']

    def run_one_cycle(self):
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

        intensity = self.objective_function.get_values(selection)
        #print intensity
        new_selection = {}
        for imember in selection:
            new_selection[imember] = None

        # Move all the fireflies (Except the most brightness)
        for i in range(len(selection)):
            imember = selection[i]
            print imember, self.population.member_str(imember)
            for j in range(len(selection)):
                jmember = selection[j]
                distance = self.population.distance(imember, jmember)
                if abs(intensity[imember]) < 1E-7:
                    intensity[imember] += 1E-7
                if abs(intensity[jmember]) < 1E-7:
                    intensity[jmember] += 1E-7
                # Minus sign because we are searching for minima
                if -intensity[imember] < -math.exp(-self.gamma * distance) * intensity[jmember]:
                    if new_selection[imember] is None:
                        new_selection[imember] = self.population.move(imember, jmember, in_place=False)
                    else:
                        self.population.move(new_selection[imember], jmember, in_place=True)

        for imember in sorted_selection:
            if new_selection[imember] is not None:
                self.population.disable(imember)
                new_member = new_selection[imember]
                self.generation[new_member] = [self.current_generation + 1]
                print imember, '   Changed   %5.2f   %s -> %s' % (intensity[imember], imember, new_member)
            else:
                print imember, '   Unchanged %5.2f' % intensity[imember]
                self.pass_to_new_generation(imember)

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
