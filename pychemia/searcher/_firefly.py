__author__ = 'Guillermo Avendano-Franco'

import math
from pychemia import log
from _searcher import Searcher


class FireFly(Searcher):

    def __init__(self, population, params, fraction_evaluated=0.8, generation_size=32, stabilization_limit=10):
        """
        Implementation fo the Firefly algorithm for global minimization
        This searcher uses a metric to compute the attractiveness and the vector displacement
        to move one firefly in the direction of another one

        :param population:
        :param params: (dict) Parameters to setup the Searcher
        :param number_generations:
        :param fraction_evaluated:
        :return:
        """
        # Mandatory objects
        self.population = population
        # Parameters
        self.gamma = None
        self.set_params(params)
        # Constrains
        self.fraction_evaluated = fraction_evaluated
        self.generation_size = generation_size
        self.stabilization_limit = stabilization_limit
        # Initializing objects
        Searcher.__init__(self, self.population, fraction_evaluated, generation_size, stabilization_limit)

    def set_params(self, params):
        assert('gamma' in params)
        assert(params['gamma'] >= 0.0)
        self.gamma = params['gamma']

    def run_one_cycle(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.get_generation_evaluated()
        sorted_selection = self.population.ids_sorted(selection)
        print 'Size of selection : ', len(selection)
        print 'Size of actives   : ', len(self.population.actives)
        print 'Size of members   : ', len(self.population.members)
        print 'Size of generation   : ', len(self.generation)
        self.print_status()

        # Automatic promotion for active members that are not evaluated
        for imember in [x for x in self.population.actives if x not in sorted_selection]:
            print imember, " Active, not evaluated, promoted "
            self.pass_to_new_generation(imember)

        intensity = self.population.get_values(selection)
        new_selection = {}
        for imember in selection:
            new_selection[imember] = None

        # Move all the fireflies (Except the most brightness)
        for i in range(len(selection)):
            imember = selection[i]
            print imember
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

