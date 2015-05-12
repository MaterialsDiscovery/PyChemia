__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta, abstractmethod
from pychemia import pcm_log
import time


class Searcher():
    __metaclass__ = ABCMeta
    """
    General class for all optimization algorithms that uses fixed and blocked
    Generations
    """

    def __init__(self, population, fraction_evaluated=0.9, generation_size=32, stabilization_limit=10):
        self.population = population
        self.generation_size = generation_size
        self.fraction_evaluated = fraction_evaluated
        self.stabilization_limit = stabilization_limit
        self.current_generation = 0
        self.generation = {}
        self.sleep_time = 2
        self.delta_change = 0.2

        self.old_actives = []
        self.old_nextgen = []
        self.pcdb = population.pcdb

    def recover(self):
        if self.pcdb is not None:
            data = self.pcdb.db.searcher_info.find_one({'_id': self.population.tag})
            if data is not None:
                self.fraction_evaluated = data['fraction_evaluated']
                self.stabilization_limit = data['stabilization_limit']
                self.generation_size = data['generation_size']
                self.set_params(data['params'])
                self.generation = {}
                self.current_generation = 0

            for entry in self.population.get_all_generations():
                if self.population.tag in entry:
                    self.generation[entry['_id']] = entry[self.population.tag]
                    if max(self.generation[entry['_id']]) > self.current_generation:
                        self.current_generation = max(self.generation[entry['_id']])

    def __str__(self):
        ret = ' Searcher Name:       %s\n' % self.searcher_name
        ret += ' Fraction evaluated:  %5.3f\n' % self.fraction_evaluated
        ret += ' Generation size:     %d\n' % self.generation_size
        ret += ' Stabilization limit: %d\n' % self.stabilization_limit
        ret += ' Parameters: %s\n' % str(self.get_params())
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

            if len(actives) == 0:
                self.population.random_population(self.generation_size)
                break
            elif len(actives) == self.generation_size:
                pcm_log.debug('Equal')
                break
            elif len(actives) > self.generation_size:
                pcm_log.debug('More %d' % len(actives))
                # Overpopulated removing some members
                if len(actives_evaluated) > self.generation_size:
                    n = len(actives_evaluated) - self.generation_size
                    pcm_log.debug('Disabling the %d worst entries' % n)
                    candidates = self.population.ids_sorted(actives_evaluated)
                    for entry_id in candidates[-n:]:
                        self.population.disable(entry_id)
                else:
                    n = len(actives) - self.generation_size
                    pcm_log.debug('Disabling %d non evaluated entries' % n)
                    for i in range(n):
                        self.population.disable(actives_no_evaluated[i])
            else:
                pcm_log.debug('Less %d' % len(actives))
                if len(evaluated) >= self.generation_size:
                    pcm_log.debug('Enabling some entries starting from the best candidates')
                    candidates = self.population.ids_sorted(evaluated)
                    index = 0
                    while len(actives) < self.generation_size and index < len(evaluated):
                        if candidates[index] not in actives:
                            pcm_log.debug('Not in actives %s' % str(candidates[index]))
                            self.population.enable(candidates[index])
                        else:
                            pcm_log.debug('In actives %s' % str(candidates[index]))
                        index += 1
                        actives = self.population.actives
                else:
                    n = self.generation_size - len(actives)
                    pcm_log.debug('Creating %d random structures' % n)
                    self.population.random_population(n)

    def get_generation(self, value=None):
        """
        Return all the elements tagged as belonging to a given generation

        :return:
        """
        if value is None:
            value = self.current_generation
        return [x for x in self.generation if value in self.generation[x]]

    def pass_to_new_generation(self, entry_id, reason=None):
        pcm_log.debug('Moving to new generation : %s' % str(entry_id))
        change = {'change': 'promoted', 'reason': reason}
        self.set_generation(entry_id, self.current_generation + 1)
        self.write_change(entry_id, change)

    def print_status(self, level='DEBUG'):
        if level=='DEBUG':
            pcm_log.debug(' %s (tag: %s)' % (self.population.name, self.population.tag))
            pcm_log.debug(' Current Generation             : %4d' % self.current_generation)
            pcm_log.debug(' Population (evaluated/total)   : %4d /%4d' % (len(self.population.evaluated),
                                                                    len(self.population.members)))
            pcm_log.debug(' Actives (evaluated/total)      : %4d /%4d' % (len(self.population.actives_evaluated),
                                                                    len(self.population.actives)))
            pcm_log.debug(' Size of Generation (this/next) : %4d /%4d' % (len(self.get_generation()),
                                                                len(self.get_generation(self.current_generation+1))))
        else:
            print '-'
            pcm_log.info(' %s (tag: %s)' % (self.population.name, self.population.tag))
            pcm_log.info(' Current Generation             : %4d' % self.current_generation)
            pcm_log.info(' Population (evaluated/total)   : %4d /%4d' % (len(self.population.evaluated),
                                                                    len(self.population.members)))
            pcm_log.info(' Actives (evaluated/total)      : %4d /%4d' % (len(self.population.actives_evaluated),
                                                                    len(self.population.actives)))
            pcm_log.info(' Size of Generation (this/next) : %4d /%4d' % (len(self.get_generation()),
                                                                len(self.get_generation(self.current_generation+1))))

        if len(self.get_generation(self.current_generation+1)) + len(self.population.actives) != self.generation_size:
            pcm_log.debug('Change in generations')
            for i in self.old_actives:
                if i not in self.population.actives:
                    pcm_log.debug('This active disappeared: %s' % str(i))
            for i in self.population.actives:
                if i not in self.old_actives:
                    pcm_log.debug('This active appeared: %s' % str(i))
            for i in self.old_nextgen:
                if i not in self.get_generation(self.current_generation+1):
                    pcm_log.debug('This nextgen disappeared: %s' % str(i))
            for i in self.get_generation(self.current_generation+1):
                if i not in self.old_nextgen:
                    pcm_log.debug('This nextgen appeared: %s' % str(i))
        self.old_actives = self.population.actives
        self.old_nextgen = self.get_generation(self.current_generation+1)

    def replace_by_changed(self, entry_id_old, reason=None):
        entry_id_new = self.population.move_random(entry_id_old, factor=self.delta_change, in_place=False, kind='change')
        change = {'change': 'modified', 'to': entry_id_new, 'reason': reason}
        pcm_log.debug('Modified  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        self.population.disable(entry_id_old)
        self.set_generation(entry_id_new, self.current_generation + 1)
        self.write_change(entry_id_old, change)

    def replace_by_other(self, entry_id_old, entry_id_new, reason=None):
        change = {'change': 'replace_by_other', 'to': entry_id_new, 'reason': reason}
        pcm_log.debug('Changed  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        self.population.disable(entry_id_old)
        self.set_generation(entry_id_new, self.current_generation + 1)
        self.write_change(entry_id_old, change)

    def replace_by_random(self, entry_id, reason=None):
        self.population.disable(entry_id)
        new_member, origin = self.population.add_random()
        pcm_log.debug('Replace by random  %s -> %s' % (str(entry_id), str(new_member)))
        self.set_generation(new_member, self.current_generation + 1)
        change = {'change': 'replace_by_random', 'to': new_member, 'reason': reason, 'origin': origin}
        self.write_change(entry_id, change)

    def save_generations(self):
        if self.pcdb is not None:
            for entry_id in sorted(self.generation):
                info = self.generation[entry_id]
                if self.pcdb.db.generations.find_one({'_id': entry_id}) is None:
                    self.pcdb.db.generations.insert({'_id': entry_id, self.population.tag: info})
                else:
                    self.pcdb.db.generations.update({'_id': entry_id}, {'$set': {self.population.tag: info}})

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
        print str(self.population)
        self.save_info()
        self.population.save_info()

        while True:
            self.print_status(level='DEBUG')

            pcm_log.debug('[%s] Enforcing the size of generation: %d' % (self.searcher_name, self.generation_size))
            self.enforce_generation_size()

            pcm_log.debug('Setting generation...')
            for entry_id in self.population.actives:
                self.set_generation(entry_id, self.current_generation)

            self.old_actives = self.population.actives
            self.old_nextgen = self.get_generation(self.current_generation+1)

            number_evaluated = len(self.population.actives_evaluated)
            while len(self.population.actives_evaluated) < self.fraction_evaluated * self.generation_size:
                if len(self.population.actives_evaluated) != number_evaluated:
                    pcm_log.debug("Population '%s' still not evaluated. Ratio %5.3f / %5.3f" % (self.population.name,
                                                                                            self.population.fraction_evaluated,
                                                                                            self.fraction_evaluated))
                    self.print_status(level='DEBUG')
                    number_evaluated = len(self.population.actives_evaluated)
                self.population.replace_failed()
                time.sleep(self.sleep_time)
            pcm_log.debug("Population '%s' evaluated. Ratio %5.3f / %5.3f" %
                      (self.population.name, self.population.fraction_evaluated, self.fraction_evaluated))

            pcm_log.debug('[%s] Removing not evaluated: %d' %
                     (self.searcher_name, len(self.population.actives_no_evaluated)))
            for entry_id in self.population.actives_no_evaluated:
                self.replace_by_random(entry_id, reason='no evaluated')
            self.print_status()

            duplicates = self.population.check_duplicates(self.population.actives_evaluated)
            pcm_log.debug('[%s] Removing duplicates: %d' % (self.searcher_name, len(duplicates)))
            for entry_id in duplicates:
                change = {'change': 'duplicate', 'to': duplicates[entry_id], 'reason': None}
                self.write_change(entry_id, change)
                self.replace_by_random(entry_id, reason='duplicate')
            self.print_status(level='INFO')

            pcm_log.debug('[%s] Running one cycle' % self.searcher_name)
            self.run_one()
            self.print_status(level='DEBUG')

            best_member = self.population.ids_sorted(self.population.actives_evaluated)[0]
            pcm_log.info('Best candidate: [%s] %s' % (best_member, self.population.str_entry(best_member)))
            pcm_log.debug('[%s] Generations %s' % (str(best_member), str(self.generation[best_member])))
            if len(self.generation[best_member]) > self.stabilization_limit:
                break

            # Increase the current generation number
            self.current_generation += 1

            # Enable all entries in the new generation
            # Disable all entries not in new generation
            pcm_log.debug('[%s] Activating new generation, disabling others' % self.searcher_name)
            for entry_id in self.generation:
                if self.current_generation in self.generation[entry_id]:
                    self.population.enable(entry_id)
                else:
                    self.population.disable(entry_id)
            self.print_status(level='DEBUG')

            self.save_generations()

        print 'Searcher ended after %d iterations' % self.current_generation
        print 'Best candidate: [%s] %s' % (best_member, self.population.str_entry(best_member))

    def write_change(self, entry_id, change):
        if self.pcdb is not None:
            change['method'] = self.searcher_name
            change['tag'] = self.population.tag
            change['generation'] = self.current_generation
            change['from'] = entry_id
            self.pcdb.db.generation_changes.insert(change)

    @property
    def searcher_name(self):
        return self.__class__.__name__

    @property
    def to_dict(self):
        return {'name': self.searcher_name,
                'fraction_evaluated': self.fraction_evaluated,
                'generation_size': self.generation_size,
                'stabilization_limit': self.stabilization_limit,
                'params': self.get_params()}

    def save_info(self):
        if self.pcdb is not None:
            data = self.pcdb.db.searcher_info.find_one({'_id': self.population.tag})
            if data is None:
                data = self.to_dict
                data['_id'] = self.population.tag
                self.pcdb.db.searcher_info.insert(data)
            else:
                self.pcdb.db.searcher_info.update({'_id': self.population.tag}, self.to_dict)

    def replace_failed(self):
        self.pcdb.replace_failed()

    def get_all_generations(self):
        return self.pcdb.db.generations.find()
