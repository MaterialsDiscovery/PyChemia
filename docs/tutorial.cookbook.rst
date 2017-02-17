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

Converting an ascii file into a POSCAR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider the following ascii file stored on ``test.ascii``:

.. code-block:: none

    None
    11.47012476778924 0.79937702290141 9.51246277071292
    -2.99939838492446 -0.12947182393907 7.79142604544631
       0.6793241103939325 -0.0865078526005124  5.2059167421975845  Mg
       7.2961243488047218  7.8694796435395364  5.3644840894866279  Mg
       2.9186811537636199  0.3979635592494910  1.4510323645408516  Ca
      -0.9269026176064651  2.5195450955593404  2.9685499969876323  Ca
      -0.8684364924472798  6.8223765997007311  5.7332541305138323  Ca
       7.7311137749974712  6.1831255116012693  0.0000000000000000  Ca

This is how pychemia can convert it into a POSCAR file:


    >>> import pychemia
    >>> st = pychemia.io.ascii.load('test.ascii')
    >>> pychemia.code.vasp.write_poscar(st, 'POSCAR.test')

The final archive calle ``POSCAR.test`` will looks like this:

.. code-block:: none

     Mg Ca
    1.0
      11.4701247677892404   0.0000000000000000   0.0000000000000000
       0.7993770229014100   9.5124627707129203   0.0000000000000000
      -2.9993983849244601  -0.1294718239390700   7.7914260454463102
     Mg Ca
     2 4
    Direct
       0.2339469912681612   1.0000000000000000   0.6681596811459407
       0.7578333625215485   0.8366521516169854   0.6885111991304717
       0.3000667985403045   0.0443708102021889   0.1862345039376849
       0.0000000000000000   0.2700535286675975   0.3810021400026761
       0.0660256684506364   0.7272193856947440   0.7358414360955946
       0.6287217253130448   0.6500025977119140   0.0000000000000000
