Global Minimimization of Lennard-Jones Clusters
-----------------------------------------------

This tutorial will guide to how search for global minima using the methods implemented on PyChemia.

Quick version
~~~~~~~~~~~~~

The shortest version of a global search using the FireFly method will look like this

    >>> from pychemia.searcher import FireFly
    >>> from pychemia.population import LJCluster
    >>> popu = LJCluster('LJ13', composition='Xe13', refine=True, direct_evaluation=True)
    >>> searcher = FireFly(popu, generation_size=16, stabilization_limit=10)
    >>> searcher.run()

For this case, you should have a mongo server running on you local machine, no SSL encryption
and no authorization with username and password.
The population will be created with Lennard-Jones clusters with 13 particles each.
Each new candidate is locally relaxed when created.
The searcher will use 16 candidates on each generation and will stop when the best candidate
survives for 10 generations.

