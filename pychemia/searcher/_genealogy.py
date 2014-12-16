__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta


class Genealogy():
    __metaclass__ = ABCMeta
    """
    Abstract class for all optimization algorithms that uses MetaHeuristics
    """

    def __init__(self, population, evaluator):
        self.population = population
        self.evaluator = evaluator
        self.current_generation = 0
        self.generation = {}
        for i in self.population.actives:
            self.generation[i] = [self.current_generation]

    def ratio_evaluated_active(self):
        """
        Return the fraction of elements considered as relaxed
        compare with the number of elements

        :return:
        """
        generation = self.get_generation()
        evaluated = self.evaluator.evaluated
        actives = self.population.actives

        number_of_evaluated_and_active = 0
        for i in generation:
            if i in evaluated and i in actives:
                number_of_evaluated_and_active += 1

        return number_of_evaluated_and_active / len(generation)

    def get_generation(self):
        """
        Return all the elements tagged as belonging to the current
        generation

        :return:
        """
        return [x for x in self.generation if self.current_generation in self.generation[x]]

    def get_generation_evaluated(self):
        """
        Return all the elements in the current generation that are also relaxed

        :return:
        """
        evaluated = self.population.evaluated
        #print 'Evaluated ', len(relaxed)
        #print 'Generation ', len(self.generation)
        #print 'Current Gen', self.current_generation
        #for i in sorted(self.generation):
        #    print i, self.generation[i]
        return [x for x in self.generation if self.current_generation in self.generation[x] and x in evaluated]

    def pass_to_new_generation(self, imember):
        self.generation[imember].append(self.current_generation + 1)

    def save_generations(self):
        wf = open('generations.dat', 'w')
        for i in sorted(self.generation):
            wf.write('%15s %10d\n' % (i, len(self.generation[i])))
        wf.close()

    def replacing(self, imember):
        self.population.disable(imember)
        new_member = self.population.add_random()
        # Replacing the element in the same generation
        ngeneration = self.generation[imember][-1]
        self.generation[new_member] = [ngeneration]
        self.generation[imember].remove(ngeneration)

    def print_status(self):
        print 'Generation', self.current_generation, sorted(self.get_generation())
        print 'Selection ', self.current_generation, sorted(self.get_generation_evaluated())
        print 'Actives   ', self.current_generation, sorted(self.population.actives)

    def __len__(self):
        return len(self.population)

