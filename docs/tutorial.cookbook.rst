PyChemia CookBook
-----------------

This is a set of small recipes for getting the 'job done' quickly
using the methods and classes implemented on PyChemia

Reading a formula and extracting the composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider that you have a formula such as 'YBa2Cu3O7' and you would like
to get the composition and the number of atoms of each species.
The recipe is very simple:

    >>> import pychemia
    >>> formula = 'YBa2Cu3O7'
    >>> comp = pychemia.Composition(formula)
    >>> comp['Ba']
    2

The Composition object acts like a python dictionary, with the particularity of returning 0 when we ask for species
non present on a given composition:

    >>> comp['Au']
    0


