__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta, abstractmethod
from pychemia import log
import time
import json


class Searcher():
    __metaclass__ = ABCMeta
    """
    Abstract class for all optimization algorithms that uses fixed and blocked
    Generations
    """

    def __init__(self, population, fraction_evaluated=0.95, generation_size=32, stabilization_limit=10):
        self.population = population
        self.generation_size = generation_size
        self.fraction_evaluated = fraction_evaluated
        self.stabilization_limit = stabilization_limit
        self.current_generation = 0
        self.generation = {}
        self.generation_changes = {}
        self.sleep_time = 2

        self.old_actives = []
        self.old_nextgen = []

    def __str__(self):
        ret = ' Searcher Name:       %s\n' % self.searcher_name
        ret += ' Fraction evaluated:  %5.3f\n' % self.fraction_evaluated
        ret += ' Generation size:     %d\n' % self.generation_size
        ret += ' Stabilization limit: %d\n\n' % self.stabilization_limit
        ret += str(self.population)
        return ret

    def enforce_generation_size(self):

        while True:
            # Get a snapshot of the current status
            # Otherwise some inconsistencies could appear if some entries become evaluated during the execution
            # of this routine
            actives = self.population.actives
            evaluated = self.population.evaluated
            actives_evaluated = self.population.actives_evaluated
            actives_no_evaluated = self.population.actives_no_evaluated

            if len(actives) == self.generation_size:
                log.debug('Equal')
                break
            elif len(actives) > self.generation_size:
                log.debug('More %d' % len(actives))
                # Overpopulated removing some members
                if len(actives_evaluated) > self.generation_size:
                    n = len(actives_evaluated) - self.generation_size
                    log.debug('Disabling the %d worst entries' % n)
                    candidates = self.population.ids_sorted(actives_evaluated)
                    for entry_id in candidates[-n:]:
                        self.population.disable(entry_id)
                else:
                    n = len(actives) - self.generation_size
                    log.debug('Disabling %d non evaluated entries' % n)
                    for i in range(n):
                        self.population.disable(actives_no_evaluated[i])
            else:
                log.debug('Less %d' % len(actives))
                if len(evaluated) >= self.generation_size:
                    log.debug('Enabling some entries starting from the best candidates')
                    candidates = self.population.ids_sorted(evaluated)
                    index = 0
                    while len(actives) < self.generation_size and index < len(evaluated):
                        if candidates[index] not in actives:
                            log.debug('Not in actives %s' % str(candidates[index]))
                            self.population.enable(candidates[index])
                        else:
                            log.debug('In actives %s' % str(candidates[index]))
                        index += 1
                        actives = self.population.actives
                else:
                    n = self.generation_size - len(actives)
                    log.debug('Creating %d random structures' % n)
                    self.population.random_population(n)

    def get_generation(self, value=None):
        """
        Return all the elements tagged as belonging to a given generation

        :return:
        """
        if value is None:
            value = self.current_generation
        return [x for x in self.generation if value in self.generation[x]]

    # @property
    # def generation_evaluated(self):
    #     """
    #     Return all the elements in the current generation that are also relaxed
    #
    #     :return:
    #     """
    #     return [x for x in self.get_generation() if self.population.is_evaluated(x)]
    #
    # @property
    # def generation_no_evaluated(self):
    #     """
    #     Return all the elements in the current generation that are also relaxed
    #
    #     :return:
    #     """
    #     return [x for x in self.get_generation() if not self.population.is_evaluated(x)]

    def pass_to_new_generation(self, entry_id, reason=None):
        log.debug('Moving to new generation : %s' % str(entry_id))
        change = {'change': 'promoted', 'reason': reason}
        self.set_generation(entry_id, self.current_generation + 1)
        self.write_change(entry_id, change)

    def print_status(self):
        print '\n'
        log.info(' GENERATION           : %4d' % self.current_generation)
        log.info(' Population           : %4d' % len(self.population.members))
        log.info(' Actives              : %4d' % len(self.population.actives))
        log.info(' Actives evaluated    : %4d' % len(self.population.actives_evaluated))
        log.info(' This Generation      : %4d' % len(self.get_generation()))
        log.info(' Next Generation      : %4d' % len(self.get_generation(self.current_generation+1)))

        if len(self.get_generation(self.current_generation+1)) + len(self.population.actives) != self.generation_size:
            log.critical('Anomalous situation')
            for i in self.old_actives:
                if i not in self.population.actives:
                    log.critical('This active disappeared: %s' % str(i))
            for i in self.population.actives:
                if i not in self.old_actives:
                    log.critical('This active appeared: %s' % str(i))
            for i in self.old_nextgen:
                if i not in self.get_generation(self.current_generation+1):
                    log.critical('This nextgen disappeared: %s' % str(i))
            for i in self.get_generation(self.current_generation+1):
                if i not in self.old_nextgen:
                    log.critical('This nextgen appeared: %s' % str(i))
        self.old_actives = self.population.actives
        self.old_nextgen = self.get_generation(self.current_generation+1)

    def replace_by_changed(self, entry_id_old, reason=None):
        entry_id_new = self.population.add_modified(entry_id_old)
        change = {'change': 'modified', 'to': entry_id_new, 'reason': reason}
        log.debug('Modified  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        self.population.disable(entry_id_old)
        self.set_generation(entry_id_new, self.current_generation + 1)
        self.write_change(entry_id_old, change)

    def replace_by_other(self, entry_id_old, entry_id_new, reason=None):
        change = {'change': 'replace_by_other', 'to': entry_id_new, 'reason': reason}
        log.debug('Changed  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        self.population.disable(entry_id_old)
        self.set_generation(entry_id_new, self.current_generation + 1)
        self.write_change(entry_id_old, change)

    def replace_by_random(self, entry_id, reason=None):
        self.population.disable(entry_id)
        new_member = self.population.add_random()
        log.debug('Replace by random  %s -> %s' % (str(entry_id), str(new_member)))
        self.set_generation(new_member, self.current_generation + 1)
        change = {'change': 'replace_by_random', 'to': new_member, 'reason': reason}
        self.write_change(entry_id, change)

    def save_generations(self):
        wf = open('generations.json', 'w')
        ret = {}
        for entry_id in sorted(self.generation):
            ret[str(entry_id)] = self.generation[entry_id]

            if self.population.pcdb.db.generations.find_one({'_id': entry_id}) is None:
                self.population.pcdb.db.generations.insert({'_id': entry_id,
                                                            self.population.tag: self.generation[entry_id]})
            else:
                self.population.pcdb.db.generations.update({'_id': entry_id},
                                                           {'$set': {self.population.tag: self.generation[entry_id]}})

        json.dump(ret, wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()

    def save_generation_changes(self):
        wf = open('generation_changes.txt', 'w')
        wf.write(str(self.generation_changes))
        wf.close()

    def set_generation(self, entry_id, value):
        if entry_id not in self.generation:
            self.generation[entry_id] = [value]
        else:
            if value not in self.generation[entry_id]:
                self.generation[entry_id].append(value)

    @abstractmethod
    def run_one(self):
        pass

    @abstractmethod
    def set_params(self, params):
        pass

    @abstractmethod
    def get_params(self):
        pass

    def run(self):
        """
        Execute the total number of cycles

        :return:
        """
        print str(self)
        self.save_info()
        self.population.save_info()

        while True:
            self.print_status()

            log.info('[%s] Enforcing the size of generation: %d' % (self.searcher_name, self.generation_size))
            self.enforce_generation_size()

            for entry_id in self.population.actives:
                self.set_generation(entry_id, self.current_generation)

            self.old_actives = self.population.actives
            self.old_nextgen = self.get_generation(self.current_generation+1)

            number_evaluated = len(self.population.actives_evaluated)
            while len(self.population.actives_evaluated) < self.fraction_evaluated * self.generation_size:
                if len(self.population.actives_evaluated) != number_evaluated:
                    log.debug("Population '%s' still not evaluated. Ratio %5.3f / %5.3f" % (self.population.name,
                                                                                            self.population.fraction_evaluated,
                                                                                            self.fraction_evaluated))
                    self.print_status()
                    number_evaluated = len(self.population.actives_evaluated)
                time.sleep(self.sleep_time)
            log.debug("Population '%s' evaluated. Ratio %5.3f / %5.3f" %
                      (self.population.name, self.population.fraction_evaluated, self.fraction_evaluated))

            log.info('[%s] Removing not evaluated: %d' %
                     (self.searcher_name, len(self.population.actives_no_evaluated)))
            for entry_id in self.population.actives_no_evaluated:
                self.replace_by_random(entry_id, reason='no evaluated')
            self.print_status()

            duplicates = self.population.check_duplicates(self.population.actives_evaluated)
            log.info('[%s] Removing duplicates: %d' % (self.searcher_name, len(duplicates)))
            for entry_id in duplicates:
                change = {'change': 'duplicate', 'to': duplicates[entry_id], 'reason': None}
                self.write_change(entry_id, change)
                self.replace_by_random(entry_id, reason='duplicate')
            self.print_status()

            log.info('[%s] Running one cycle' % self.searcher_name)
            self.run_one()
            self.print_status()

            best_member = self.population.ids_sorted(self.population.actives_evaluated)[0]
            log.info('[%s] Best member is %s' % (self.searcher_name, str(best_member)))
            log.debug('[%s] Generations %s' % (str(best_member), str(self.generation[best_member])))
            if len(self.generation[best_member]) > self.stabilization_limit:
                break

            # Increase the current generation number
            self.current_generation += 1

            # Enable all entries in the new generation
            # Disable all entries not in new generation
            log.info('[%s] Activating new generation, disabling others' % self.searcher_name)
            for entry_id in self.generation:
                if self.current_generation in self.generation[entry_id]:
                    self.population.enable(entry_id)
                else:
                    self.population.disable(entry_id)
            self.print_status()

            self.save_generations()
            self.save_generation_changes()

    def write_change(self, entry_id, change):
        change['method'] = self.searcher_name
        change['tag'] = self.population.tag
        change['generation'] = self.current_generation
        change['from'] = entry_id
        self.population.pcdb.db.generation_changes.insert(change)
        if entry_id not in self.generation_changes:
            self.generation_changes[entry_id] = [change]
        else:
            self.generation_changes[entry_id].append(change)

    @property
    def searcher_name(self):
        return self.__class__.__name__

    def to_dict(self):
        return {'name': self.searcher_name,
                'fraction_evaluated': self.fraction_evaluated,
                'generation_size': self.generation_size,
                'stabilization_limit': self.stabilization_limit,
                'params': self.get_params()}

    def save_info(self):
        self.population.pcdb.db.searcher_info.insert(self.to_dict())
