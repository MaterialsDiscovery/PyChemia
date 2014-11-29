__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta


class MetaHeuristics():
    __metaclass__ = ABCMeta
    """
    Abstract class for all optimization algorithms that uses MetaHeuristics
    """

    def __init__(self, population, relax_population):
        self.population = population
        self.relax_population = relax_population
        self.current_generation = 0
        self.generation = {}
        for i in self.population.all_entries:
            self.generation[i] = [self.current_generation]

    def ratio_relaxed_active(self):
        """
        Return the fraction of elements considered as relaxed
        compare with the number of elements

        :return:
        """
        generation = self.get_generation()
        relaxed = self.relax_population.relaxed
        actives = self.population.actives

        number_of_relaxed_and_active = 0
        for i in generation:
            if i in relaxed and i in actives:
                number_of_relaxed_and_active += 1

        return number_of_relaxed_and_active / len(generation)

    def get_generation(self):
        """
        Return all the elements tagged as belonging to the current
        generation

        :return:
        """
        return [x for x in self.generation if self.current_generation in self.generation[x]]

    def get_generation_relaxed(self):
        """
        Return all the elements in the current generation that are also relaxed

        :return:
        """
        relaxed = self.relax_population.relaxed
        return [x for x in self.generation if self.current_generation in self.generation[x] and x in relaxed]

    def pass_to_new_generation(self, imember):
        self.generation[imember].append(self.current_generation + 1)
