
import time
from abc import ABCMeta, abstractmethod
from pychemia import pcm_log, HAS_PYMONGO

if HAS_PYMONGO:
    import bson


class Searcher:
    __metaclass__ = ABCMeta
    """
    General class for all optimization algorithms that uses fixed and blocked
    Generations
    """

    def __init__(self, population, generation_size=32, stabilization_limit=10,
                 target_value=None, searcher_id=None):
        self.population = population
        self.generation_size = generation_size
        self.stabilization_limit = stabilization_limit
        self.target_value = target_value
        self.current_generation = 0
        self.generation = {}
        self.sleep_time = 2
        self.delta_change = 0.2
        self.lineage = {}
        self.lineage_inv = {}

        self.old_actives = []
        self.old_nextgen = []
        self.pcdb = population.pcdb
        if searcher_id is None:
            self.searcher_id = self.population.tag
        else:
            self.searcher_id = searcher_id

    def recover(self, changedb=False):
        if self.pcdb is not None:
            data = self.pcdb.db.searcher_info.find_one({'_id': self.searcher_id})
            if data is not None:
                self.stabilization_limit = data['stabilization_limit']
                self.generation_size = data['generation_size']
                self.set_params(data['params'])
                self.generation = {}
                self.current_generation = 0

            for entry in self.get_all_generations():
                if self.searcher_id in entry:
                    self.generation[entry['_id']] = entry[self.searcher_id]
                    if max(self.generation[entry['_id']]) > self.current_generation:
                        self.current_generation = max(self.generation[entry['_id']])

            lineage_data = self.pcdb.db.lineage.find_one({'_id': self.searcher_id})
            if lineage_data is not None:
                lineage = lineage_data['lineage']
                lineage_inv = lineage_data['lineage_inv']

                changed = False
                for i in lineage.keys():
                    for j in lineage[i]:
                        if j == 'None':
                            lineage[i].remove(j)
                            changed = True
                if changed and changedb:
                    print('Lineage was changed (Some entry was null)')
                    self.pcdb.db.lineage.update_one({'_id': self.searcher_id}, {'$set': {'lineage': lineage}})

                sizes = [len(lineage[x]) for x in lineage]
                assert min(sizes) == max(sizes)

                for i in lineage.keys():
                    if HAS_PYMONGO:
                        self.lineage[i] = [bson.ObjectId(x) for x in lineage[i]]
                    else:
                        self.lineage[i] = [x for x in lineage[i]]

                for i in lineage_inv.keys():
                    if HAS_PYMONGO:
                        self.lineage_inv[bson.ObjectId(i)] = lineage_inv[i]
                    else:
                        self.lineage_inv[i] = lineage_inv[i]

                self.correct_extras(changedb)

    def correct_extras(self, changedb=False):
        for entry_id in self.get_generation():
            if entry_id not in self.lineage_inv:
                print('Disabling one entry not in lineage_inv', entry_id)
                self.population.disable(entry_id)
                print(self.generation.pop(entry_id))
                if changedb:
                    self.pcdb.db.generations.remove({'_id': entry_id})
            else:
                slot = self.lineage_inv[entry_id]
                if self.lineage[slot][-1] != entry_id:
                    print('Disabling one entry not in lineage[slot][-1]', entry_id)
                    self.population.disable(entry_id)
                    print(self.generation.pop(entry_id))
                    if changedb:
                        self.pcdb.db.generations.remove({'_id': entry_id})

        if self.current_generation > 0:
            for slot in range(self.generation_size):
                entry_id = self.lineage[str(slot)][-1]
                if entry_id not in self.population.actives:
                    print('Activating from lineage', entry_id)
                    self.population.enable(entry_id)
                if entry_id not in self.get_generation():
                    self.set_generation(entry_id, self.current_generation)
            actives = self.population.actives
            for entry_id in actives:
                if entry_id not in self.get_generation():
                    print('Disabling ', entry_id)
                    self.population.disable(entry_id)
            for entry_id in self.get_generation():
                if entry_id not in self.population.actives:
                    print('Enabling', entry_id)
                    self.population.enable(entry_id)
            candidates_per_generation = [len(self.get_generation(i)) for i in range(self.current_generation + 1)]
            pcm_log.info('Candidates per generation: %s' % candidates_per_generation)
            pcm_log.info('Current generation: %d Candidates: %d' % (self.current_generation,
                                                                    len(self.get_generation())))
            print('CANDIDATES:',candidates_per_generation)
            assert len(self.get_generation()) == self.generation_size
            assert min(candidates_per_generation) == max(candidates_per_generation)

    def __str__(self):
        ret = ' Searcher Name:       %s\n' % self.searcher_name
        ret += ' Generation size:     %d\n' % self.generation_size
        ret += ' Stabilization limit: %d\n' % self.stabilization_limit
        ret += ' Current Generation:  %d\n' % self.current_generation
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
                self.correct_extras(changedb=True)
                break
            elif len(actives) > self.generation_size:
                pcm_log.debug('More %d' % len(actives))
                self.recover(changedb=True)
                self.correct_extras(changedb=True)
            else:
                pcm_log.debug('Less %d' % len(actives))
                sizes = [len(self.lineage[x]) for x in self.lineage]
                assert min(sizes) == max(sizes)
                for i in self.lineage:
                    entry_id = self.lineage[i][-1]
                    print('Activating', i, entry_id)
                    self.population.enable(entry_id)
                assert len(self.population.actives) == self.generation_size
                assert len(self.get_generation()) == self.generation_size
                sizes = [len(self.get_generation(i)) for i in range(self.current_generation)]
                assert min(sizes) == max(sizes)

                # if len(evaluated) >= self.generation_size:
                #     pcm_log.debug('Enabling some entries starting from the best candidates')
                #     candidates = self.population.ids_sorted(evaluated)
                #     index = 0
                #     while len(actives) < self.generation_size and index < len(evaluated):
                #         if candidates[index] not in actives:
                #             pcm_log.debug('Not in actives %s' % str(candidates[index]))
                #             self.population.enable(candidates[index])
                #         else:
                #             pcm_log.debug('In actives %s' % str(candidates[index]))
                #         index += 1
                #         actives = self.population.actives
                # else:
                #     n = self.generation_size - len(actives)
                #     pcm_log.debug('Creating %d random structures' % n)
                #     self.population.random_population(n)

        for i in self.population.actives:
            self.set_generation(i, self.current_generation)

    def get_generation(self, generation_number=None):
        """
        Return all the elements tagged as belonging to a given generation

        :param generation_number: Number of generation to get
        :return: (list) Identifiers of all candidates for a given generation number
        :rtype: list
        """
        if generation_number is None:
            generation_number = self.current_generation
        from_gens = [x for x in self.generation if generation_number in self.generation[x]]
        ret = []
        for i in from_gens: 
            if i in self.lineage_inv:
                ret.append(i)
        return ret

    def print_status(self):
        pcm_log.info(' %s (tag: %s)' % (self.population.name, self.population.tag))
        pcm_log.info(' Current Generation             : %4d' % self.current_generation)
        pcm_log.info(' Population (evaluated/total)   : %4d /%4d' % (len(self.population.evaluated),
                                                                     len(self.population.members)))
        pcm_log.info(' Actives (evaluated/total)      : %4d /%4d' % (len(self.population.actives_evaluated),
                                                                     len(self.population.actives)))
        pcm_log.info(' Size of Generation (this/next) : %4d /%4d\n' % (len(self.get_generation()),
                                                                       len(self.get_generation(
                                                                         self.current_generation + 1))))

        if len(self.get_generation(self.current_generation + 1)) + len(self.population.actives) != self.generation_size:
            pcm_log.debug('Change in generations')
            for i in self.old_actives:
                if i not in self.population.actives:
                    pcm_log.debug('This active disappeared: %s' % str(i))
            for i in self.population.actives:
                if i not in self.old_actives:
                    pcm_log.debug('This active appeared: %s' % str(i))
            for i in self.old_nextgen:
                if i not in self.get_generation(self.current_generation + 1):
                    pcm_log.debug('This nextgen disappeared: %s' % str(i))
            for i in self.get_generation(self.current_generation + 1):
                if i not in self.old_nextgen:
                    pcm_log.debug('This nextgen appeared: %s' % str(i))
        self.old_actives = self.population.actives
        self.old_nextgen = self.get_generation(self.current_generation + 1)

    def advance(self, father, son, change):
        self.write_change(father, change)
        if father not in self.lineage_inv:
            print('%s not in current lineage' % father)
            print('Lineages %s' % self.lineage_inv.keys())
            raise ValueError('Father not in lineage')
        slot = self.lineage_inv[father]
        self.lineage[slot].append(None)
        self.lineage[slot][self.current_generation + 1] = son
        self.lineage_inv[son] = slot
        self.set_generation(son, self.current_generation + 1)

    def pass_to_new_generation(self, entry_id, reason=None):
        # pcm_log.debug('Moving to new generation : %s' % str(entry_id))
        change = {'change': 'promoted', 'reason': reason}
        self.advance(entry_id, entry_id, change)

    def replace_by_changed(self, entry_id_old, reason=None):
        entry_id_new = self.population.move_random(entry_id_old, factor=self.delta_change, in_place=False,
                                                   kind='change')
        change = {'change': 'modified', 'to': entry_id_new, 'reason': reason}
        pcm_log.debug('Modified  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        self.population.disable(entry_id_old)
        self.advance(entry_id_old, entry_id_new, change)

    def replace_by_other(self, entry_id_old, entry_id_new, reason=None):
        change = {'change': 'replace_by_other', 'to': entry_id_new, 'reason': reason}
        # pcm_log.debug('Changed  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        self.population.disable(entry_id_old)
        self.advance(entry_id_old, entry_id_new, change)

    def replace_by_random(self, entry_id_old, reason=None):
        self.population.disable(entry_id_old)
        entry_id_new, origin = self.population.add_random()
        pcm_log.debug('Replace by random  %s -> %s' % (str(entry_id_old), str(entry_id_new)))
        change = {'change': 'replace_by_random', 'to': entry_id_new, 'reason': reason, 'origin': origin}
        self.advance(entry_id_old, entry_id_new, change)

    def save_generations(self):
        if self.pcdb is not None:
            for entry_id in sorted(self.generation):
                info = self.generation[entry_id]
                if self.pcdb.db.generations.find_one({'_id': entry_id}) is None:
                    self.pcdb.db.generations.insert_one({'_id': entry_id, self.searcher_id: info})
                else:
                    self.pcdb.db.generations.update_one({'_id': entry_id}, {'$set': {self.searcher_id: info}})

            lineage = {}
            for i in self.lineage.keys():
                lineage[i] = [str(x) for x in self.lineage[i]]

            lineage_inv = {}
            for i in self.lineage_inv.keys():
                lineage_inv[str(i)] = self.lineage_inv[i]

            if self.pcdb.db.lineage.find_one({'_id': self.searcher_id}) is None:
                self.pcdb.db.lineage.insert_one({'_id': self.searcher_id,
                                             'lineage_inv': lineage_inv,
                                             'lineage': lineage})
            else:
                self.pcdb.db.lineage.update_one({'_id': self.searcher_id}, {'$set': {'lineage': lineage,
                                                                                 'lineage_inv': lineage_inv}})

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

    def update_lineages(self):

        assert len(self.population.actives) == self.generation_size

        # Complete lineage dictionary with generation_size entries
        # Each entry on the dictionary is a list of current_generation+1 values
        for i in range(self.generation_size):
            if str(i) not in self.lineage:
                self.lineage[str(i)] = [None]
            while len(self.lineage[str(i)]) < self.current_generation + 1:
                self.lineage[str(i)].append(None)

        # Complete lineage for all actives with known lineage_inv
        for entry_id in self.population.actives:
            if entry_id in self.lineage_inv:
                slot = self.lineage_inv[entry_id]
                self.lineage[slot][self.current_generation] = entry_id

        # Find empty slots for actives not in lineage_inv
        for entry_id in self.population.actives:
            if entry_id not in self.lineage_inv:
                for i in range(self.generation_size):
                    # Search for the first empty slot and associate it to the entry_id
                    if self.lineage[str(i)][self.current_generation] is None:
                        self.lineage[str(i)][self.current_generation] = entry_id
                        self.lineage_inv[entry_id] = str(i)
                        break

        # Check consistency between lineage and lineage_inv
        for i in range(self.generation_size):
            assert str(i) in self.lineage
            assert self.lineage[str(i)][-1] in self.lineage_inv
            assert self.lineage_inv[self.lineage[str(i)][-1]] == str(i)

        # Setting Generation list
        pcm_log.debug('Setting generation...')
        for entry_id in self.population.actives:
            self.set_generation(entry_id, self.current_generation)

    def run(self):
        """
        Execute the total number of cycles

        :return:
        """

        print("Population Information")
        print("======================")
        print(self.population)
        print("Searcher Information")
        print("====================")
        print(self)
        self.save_info()
        self.population.save_info()
        best_member = ''
        best_recorded = None
        survival_for_best = 0

        while True:
            print('\nGENERATION: %d' % self.current_generation)
            if len(self.population) > 0:
                self.print_status()

            pcm_log.info('[%s] Ensuring the number of active candidates is %d' % (self.searcher_name,
                                                                                  self.generation_size))
            self.enforce_generation_size()
            self.update_lineages()

            self.old_actives = self.population.actives
            self.old_nextgen = self.get_generation(self.current_generation + 1)

            self.print_status()

            number_evaluated = len(self.population.actives_evaluated)
            while self.population.fraction_evaluated < 1.0:
                if len(self.population.actives_evaluated) != number_evaluated:
                    msg = "Population '%s' still not evaluated. %4.0f %%"
                    pcm_log.debug(msg % (self.population.name, 100 * self.population.fraction_evaluated))
                    self.print_status()
                    number_evaluated = len(self.population.actives_evaluated)
                self.population.replace_failed()
                time.sleep(self.sleep_time)
            pcm_log.info("Population '%s' evaluated. %4.0f %%" % (self.population.name,
                                                                  100 * self.population.fraction_evaluated))

            best_member = self.population.best_candidate
            self.population.refine_progressive(best_member)

            pcm_log.info('Current best candidate: [%s]:\n%s' % (best_member, self.population.str_entry(best_member)))
            if best_member in self.get_generation():
                print('This candidate have survived for %d generations' % len(self.generation[best_member]))
                if len(self.generation[best_member]) >= self.stabilization_limit:
                    self.save_generations()
                    break
            else:
                pcm_log.debug('Best candidate %s is not in the current generation' % best_member)
                # pcm_log.debug('Slot: %s' % self.lineage_inv[best_member])
                # pcm_log.debug('Lineage: %s' % self.lineage[self.lineage_inv[best_member]])
                if best_member != best_recorded:
                    survival_for_best = 0
                    best_recorded = best_member
                else:
                    survival_for_best += 1

                if survival_for_best >= self.stabilization_limit:
                    print('This candidate have survived for %d generations' % survival_for_best)
                    self.save_generations()
                    break

            if self.target_value is not None:
                if self.population.value(best_member) <= self.target_value:
                    print('Target value achieved: target=%9.3f best=%9.3f' % (self.population.value(best_member),
                                                                              self.target_value))
                    self.save_generations()
                    break
                else:
                    print('Best value = %7.3f     target value = %7.3f' % (self.population.value(best_member),
                                                                           self.target_value))

            pcm_log.debug('[%s] Removing not evaluated: %d' %
                          (self.searcher_name, len(self.population.actives_no_evaluated)))
            for entry_id in self.population.actives_no_evaluated:
                self.replace_by_random(entry_id, reason='no evaluated')
            self.print_status()

            duplicates = self.population.get_duplicates(self.population.ids_sorted(self.population.actives_evaluated))
            print(duplicates)
            for entry_id in duplicates:
                change = {'change': 'duplicate', 'to': duplicates[entry_id], 'reason': None}
                self.write_change(entry_id, change)
                self.replace_by_random(entry_id, reason='duplicate')
            pcm_log.info(' Duplicates identified and disabled: %d' % len(duplicates))
            self.print_status()

            pcm_log.info(' Running one cycle for %s with %d candidates' % (self.searcher_name,
                                                                           len(self.actives_in_generation)))
            self.run_one()
            self.update_generation()

        print('Searcher ended after %d iterations' % self.current_generation)
        print('Best candidate: [%s]:\n%s' % (best_member, self.population.str_entry(best_member)))

    def write_change(self, entry_id, change):
        if self.pcdb is not None:
            change['method'] = self.searcher_name
            change['tag'] = self.searcher_name
            change['generation'] = self.current_generation
            change['from'] = entry_id
            self.pcdb.db.generation_changes.insert_one(change)

    def update_generation(self):
        # Increase the current generation number
        assert len(self.get_generation()) == self.generation_size
        self.current_generation += 1

        # Enable all entries in the new generation
        # Disable all entries not in new generation
        pcm_log.debug('[%s] Activating new generation, disabling others' % self.searcher_name)
        for entry_id in self.generation:
            if self.current_generation in self.generation[entry_id]:
                self.population.enable(entry_id)
            else:
                self.population.disable(entry_id)
        self.save_generations()

    @property
    def searcher_name(self):
        return self.__class__.__name__

    @property
    def to_dict(self):
        return {'name': self.searcher_name,
                'generation_size': self.generation_size,
                'stabilization_limit': self.stabilization_limit,
                'params': self.get_params()}

    def save_info(self):
        if self.pcdb is not None:
            data = self.pcdb.db.searcher_info.find_one({'_id': self.searcher_id})
            if data is None:
                data = self.to_dict
                data['_id'] = self.searcher_id
                self.pcdb.db.searcher_info.insert_one(data)
            else:
                self.pcdb.db.searcher_info.update_one({'_id': self.searcher_id}, self.to_dict)

    def clean(self):
        if self.pcdb is not None:
            self.pcdb.db.searcher_info.drop()
            self.pcdb.db.generations.drop()
            self.pcdb.db.generation_changes.drop()
            self.pcdb.db.lineage.drop()

    def replace_failed(self):
        self.pcdb.replace_failed()

    def get_all_generations(self, generation_number=None):
        if generation_number is None:
            return self.pcdb.db.generations.find()
        else:
            return self.pcdb.db.generations.find({self.population.tag: generation_number})

    @property
    def actives_in_generation(self):
        return [x for x in self.population.actives if self.population.is_evaluated(x) and x in self.get_generation()]
