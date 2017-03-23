import json
import os

import bson
import networkx as nx
from networkx.readwrite import json_graph


class NetworkAnalysis:
    def __init__(self, population, tolerance=0.1):
        self.population = population
        self.tolerance = tolerance
        self.pcdb = self.population.pcdb

    def get_network(self):

        # print 'Population evaluated', len(self.population.evaluated)
        lista = list(self.population.ids_sorted(self.population.evaluated))

        graph = nx.Graph()

        for change in self.pcdb.db.generation_changes.find({'change': 'replace_by_other'}):
            graph.add_edge(str(change['from']), str(change['to']))

        for change in self.pcdb.db.generation_changes.find({'change': 'duplicate'}):
            graph.add_edge(str(change['from']), str(change['to']))

        # print 'Number of nodes:', graph.number_of_nodes()
        # print 'Number of edges:', graph.number_of_edges()

        for n in graph:
            graph.node[n]['name'] = n

        for n in graph:
            entry_id = bson.ObjectId(n)
            if 'spacegroup' in self.population.get_entry(bson.ObjectId(n))['properties']:
                graph.node[n]['group'] = str(self.population.get_entry(bson.ObjectId(n))['properties']['spacegroup'])
                graph.node[n]['name'] = graph.node[n]['group']
            else:
                graph.node[n]['group'] = str(20 * lista.index(entry_id) / len(lista))
                graph.node[n]['name'] = str(lista.index(entry_id))

        # print 'Computing distances...'
        for n in graph.edges():
            graph[n[0]][n[1]]['value'] = int(10 * self.population.distance(bson.ObjectId(n[0]), bson.ObjectId(n[1])))
        # print 'done'

        d = json_graph.node_link_data(graph)
        json.dump(d, open('Network_%s.json' % self.population.name, 'w'), sort_keys=True, indent=4,
                  separators=(',', ': '))

        if not os.path.isfile('NetworkBasins_%s.json' % self.population.name):
            dupes_dict, dupes_list = self.population.dict_duplicates(self.population.evaluated, fast=True)
            new_dupes_dict = {}
            for i in dupes_dict:
                new_dupes_dict[str(i)] = []
                for j in dupes_dict[i]:
                    new_dupes_dict[str(i)].append(str(j))

            wf = open('NetworkBasins_%s.json' % self.population.name, 'w')
            json.dump(new_dupes_dict, wf, indent=2, separators=(',', ': '))
            wf.close()
        else:
            rf = open('NetworkBasins_%s.json' % self.population.name, 'r')
            new_dupes_dict = json.load(rf)

        # print 'Number of non-duplicates', len(new_dupes_dict)

        # print 'Reverting dictionary...'
        graph_basins = nx.Graph()
        tabla_reversa = {}
        for i in lista:
            if str(i) in new_dupes_dict:
                for j in new_dupes_dict[str(i)]:
                    # print lista.index(bson.ObjectId(i)),' : ', lista.index(bson.ObjectId(j))
                    if j not in tabla_reversa:
                        tabla_reversa[j] = str(i)
        # print 'done'

        for i in graph.edges_iter():
            if i[0] in tabla_reversa and i[1] in tabla_reversa:
                graph_basins.add_edge(tabla_reversa[i[0]], tabla_reversa[i[1]])

        # print 'Number of nodes:', graph_basins.number_of_nodes()
        # print 'Number of edges:', graph_basins.number_of_edges()

        for n in graph_basins:
            entry_id = bson.ObjectId(n)
            if 'spacegroup' in self.population.get_entry(entry_id)['properties']:
                spacegroup = str(self.population.get_entry(entry_id)['properties']['spacegroup'])
                graph_basins.node[n]['group'] = spacegroup
                graph_basins.node[n]['name'] = spacegroup
            else:
                graph_basins.node[n]['group'] = str(20 * lista.index(entry_id) / len(lista))
                graph_basins.node[n]['name'] = str(lista.index(entry_id))

        for n in graph_basins.edges():
            graph_basins[n[0]][n[1]]['value'] = int(
                10 * self.population.distance(bson.ObjectId(n[0]), bson.ObjectId(n[1])))

        d2 = json_graph.node_link_data(graph_basins)
        json.dump(d2, open('Network_%s_basins.json' % self.population.name, 'w'), sort_keys=True, indent=4,
                  separators=(',', ': '))

        best = lista[0]
        line = 'Population: %10s Evaluated: %4d   N:%5d E:%5d    BN:%4d BE:%4d  %7.2f  GG: %8s %8s'
        print(line % (self.population.name,
                      len(self.population),
                      graph.number_of_nodes(),
                      graph.number_of_nodes(),
                      graph_basins.number_of_nodes(),
                      graph_basins.number_of_edges(),
                      float(graph_basins.number_of_edges()) / (1.0 + graph_basins.number_of_nodes()),
                      str(best) in graph.nodes(),
                      str(best) in graph_basins.nodes()))
