__author__ = 'Guillermo Avendano-Franco'

import math
from _metaheuristics import MetaHeuristics
import time

class FireFly(MetaHeuristics):

    def __init__(self, population, objective_function, relax_population, metric, gamma, number_generations,
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
        self.metric = metric
        self.gamma = self.set_gamma(gamma)
        self.number_generations = number_generations
        self.ratio_relax = ratio_relax
        self.objective_function.initialize(self.population)
        self.relax_population.initialize(self.population)
        self.metric.initialize(self.population)
        MetaHeuristics.__init__(self, self.population, self.relax_population)

    def __len__(self):
        return len(self.population)

    def set_gamma(self, gamma):
        assert(hmcr >= 0.0)
        self.gamma = gamma

    def run_one_cycle(self):

        # Get a static selection of the values in the generation that are relaxed
        selection = self.get_generation_relaxed()
        intensity = self.objective_function.get_values(selection)
        new_selection = {}
        for imember in selection:
            new_selection[imember] = None

        # Move all the fireflies (Except the most brightness)
        for i in range(len(selection)):
            imember = selection(i)
            for j in range(len(selection)):
                jmember = selection(j)
                distance = self.metric.distance(imember, jmember)
                if intensity[imember] < math.exp(-self.gamma * distance) * intensity[jmember]:
                    if new_selection[imember] is None:
                        new_selection[imember] = self.population.move(imember, jmember)
                    else:
                        self.population.move(new_selection[imember], jmember, in_place=True)

        for imember in selection:
            if new_selection[imember] is None:

            else:
                self.population.disable(imember)
                self.population.move(imember, jmember, in_place=new_selection[imember])


        # Increase the current generation number
        self.current_generation += 1

    def run_all_cycles(self):

        if not self.relax_population.is_running:
            self.relax_population.run()

        for i in range(self.number_generations):

            self.population.check_duplicates()

            if self.population.ratio_relaxed_active(self) > self.ratio_relax:
                self.run_one_cycle()
            else:
                time.sleep(60)
