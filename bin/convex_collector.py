import json
import argparse
import pychemia
import logging
import numpy as np

if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)

    description = """Collect structures from several databases for plotting convex hulls"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t', '--host',
                        default='localhost', metavar='server', type=str,
                        help='Hostname or address (default: localhost)')
    parser.add_argument('-o', '--port',
                        default=27017, metavar='port', type=int,
                        help='MongoDB port (default: 27017)')
    parser.add_argument('-u', '--user',
                        default=None, metavar='username', type=str,
                        help='Username (default: None)')
    parser.add_argument('-p', '--passwd',
                        default=None, metavar='password', type=str,
                        help='Password (default: None)')
    parser.add_argument('-d', '--dbname',
                        default=None, metavar='dbname', type=str, nargs='+',
                        help='PyChemia Database name (default: None)')
    parser.add_argument('-l', '--left',
                        default=None, metavar='specie', type=str,
                        help='Specie that will be on the left side')
    parser.add_argument('-r', '--right',
                        default=None, metavar='specie', type=str,
                        help='Specie that will be on the right side')
    parser.add_argument('--replicaset',
                        default=None, metavar='name', type=str,
                        help='ReplicaSet  (default: None)')
    parser.add_argument('--ssl', action='store_true',
                        help='Use SSL to connect to MongoDB  (default: No)')

    args = parser.parse_args()
    if args.dbname is None:
        parser.print_help()
        exit(1)

    print(args)
    db_settings = {'host': args.host, 'port': args.port, 'ssl': args.ssl, 'replicaset': args.replicaset}
    if args.user is not None:
        if args.passwd is None:
            raise ValueError('Password is mandatory if user is entered')
        db_settings['user'] = args.user
        db_settings['passwd'] = args.passwd

    convex_params = {'left': args.left, 'right': args.right}
    specie_right = args.right
    specie_left = args.left

    ret = []

    for idb in args.dbname:
        db_settings['name'] = idb
        pcdb = pychemia.db.get_database(db_settings)

        for entry in pcdb.entries.find():
            entry_id = entry['_id']
            st = pcdb.get_structure(entry_id)
            formula = st.formula
            symmetry = pychemia.crystal.CrystalSymmetry(st)
            space_group = symmetry.number(symprec=1e-1)

            if specie_left not in st.composition:
                ratio = 1.0
            elif specie_right not in st.composition:
                ratio = 0.0
            else:
                ratio = float(st.composition[specie_right]) / st.natom
                print("ratio: %f - right=%d natom=%d %s" % (ratio,
                                                            float(st.composition[specie_right]),
                                                            st.natom,
                                                            st.formula))

            if 'energy_pa' in entry['properties']:

                energy_pa = entry['properties']['energy_pa']

                ret.append({'formula': formula,
                            'spcgrp': space_group,
                            'energy_pa': energy_pa,
                            'ratio': ratio,
                            'natom': st.natom,
                            'density': st.density})
            else:
                print('ERROR: DB=%s ENTRY=%s No value of energy' % (idb, entry_id))

    sorter = []
    for i in ret:
        sorter.append(1E5 * i['ratio'] + i['energy_pa'])

    npsorter = np.array(sorter)
    srt = npsorter.argsort()

    newdata = []
    for i in srt:
        newdata.append(ret[i])

    wf = open('convex.json', 'w')
    json.dump(newdata, wf, sort_keys=True, indent=4, separators=(',', ': '))
    wf.close()
