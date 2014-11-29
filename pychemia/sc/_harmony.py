__author__ = 'Guillermo Avendano-Franco'

import random
from _metaheuristics import MetaHeuristics
import time


class HarmonySearch(MetaHeuristics):

    def __init__(self, population, objective_function, relax_population, hmcr, par, number_generations,
                 ratio_relax=0.8):

        self.population = population
        self.objective_function = objective_function
        self.relax_population = relax_population
        self.hmcr = self.set_hmcr(hmcr)
        self.par = self.set_par(par)
        self.number_generations = number_generations
        self.ratio_relax = ratio_relax
        self.objective_function.initialize(self.population)
        self.relax_population.initialize(self.population)
        MetaHeuristics.__init__(self, self.population, self.relax_population)

    def __len__(self):
        return len(self.population)

    def set_hmcr(self, hmcr):
        assert(1.0 >= hmcr >= 0.0)
        self.hmcr = hmcr

    def set_par(self, par):
        assert(1.0 >= par >= 0.0)
        self.par = par

    @property
    def harmony_memory_considering_rate(self):
        return self.hmcr

    def run_one_cycle(self):

        # Get a static selection of the values in the generation that are relaxed
        selection = self.get_generation_relaxed()
        for imember in selection:

            if imember in self.objective_function.top(selection):
                self.pass_to_new_generation(imember)
            if imember in self.objective_function.tail(selection):
                self.population.disable(imember)
                new_member = self.population.add_random()
                self.generation[new_member] = [self.current_generation + 1]
            else:
                rnd = random.random()
                if rnd < self.hmcr:
                    rnd = random.random()
                    if rnd < self.par:
                        new_member = self.population.add_modified(imember)
                        self.generation[new_member] = [self.current_generation + 1]
                    else:
                        self.pass_to_new_generation(imember)
                else:
                    self.population.disable(imember)
                    new_member = self.population.add_random()
                    self.generation[new_member] = [self.current_generation + 1]

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
