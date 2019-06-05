Euler Angles in PyChemia
------------------------

Short Version
~~~~~~~~~~~~~

You can use the following routines to obtain the :math:`k(k-1)/2` Generalized Euler angles
from a orhogonal matrix of dimension k and for building the orthogonal matrix from a set
of angles::

    >>> angles_list = pychemia.utils.mathematics.gea_all_angles(ortho_matrix)


    >>> ortho_matrix = gea_orthogonal_from_angles(angles_list)

Remember that the list of angles matters. The SO(k) group is non-abelian.

The algorithm is based on the paper::

    Generalization of Euler Angles to N-Dimensional Orthogonal Matrices
    David K. Hoffman, Richard C. Raffenetti, and Klaus Ruedenberg
    Journal of Mathematical Physics 13, 528 (1972)
    doi: 10.1063/1.1666011


Long Version
~~~~~~~~~~~~

Rotation on a two-dimensional plane can be described with one angle.
For a three dimensional space 3 angles are needed. You can define angles
around each axis and general rotation matrices as a product of three
single axis rotation matrices.

In general, for k dimensions, the number of angles is :math:`k(k-1)/2`.
One angle for each plane that you can get from any pair among the
k vectors defining the space.

Any orthogonal matrix (with determinant equal to +1) represent a rotation matrix
for the space that its dimension. Sometimes could be necessary to obtain the set of
independent angles that the orthogonal matrix represents and having a way
of regaining the original orthogonal matrix from a given set of angles.

One the reasons for moving between one orthogonal matrix and its angles is
reducing an orthogonal matrix to its minimal independent parameters.
The paper intitled: "Generalization of Euler Angles to N-Dimensional Orthogonal Matrices"
provides an effective way to get the so called Euler angles for matrices of arbitrary
dimension and reconstitute the matrix from a given set of angles.

This algorithm is implemented on PyChemia.
As an example, lets build first a 7-dimensional orthogonal matrix by a QR decomposition::

    In [1]: import pychemia

    In [2]: import numpy as np

    In [3]: np.set_printoptions(linewidth=200, suppress=True, precision=5)

    In [4]: ortho = pychemia.utils.mathematics.gram_smith_qr(7)

    In [5]:ortho
    Out[5]:

    array([[-0.23503793, -0.6039233 ,  0.31346146,  0.56383925,  0.18441349, -0.35292927,  0.07275729],
           [-0.31380685, -0.21108679, -0.25878856, -0.40710624, -0.23657101, -0.54173578, -0.52423003],
           [-0.48237691,  0.1056854 , -0.2376831 ,  0.11195934,  0.63817833,  0.3621268 , -0.38562619],
           [-0.44643339, -0.33091574,  0.37743259, -0.27151849, -0.39852265,  0.56178535,  0.02431587],
           [-0.46408875,  0.67693995,  0.31398053,  0.30715275, -0.28175741, -0.23147005, -0.02194835],
           [-0.43149085,  0.00696596, -0.2664925 , -0.31313612,  0.19768295, -0.20337001,  0.75117024],
           [-0.11282486, -0.10838891, -0.68280287,  0.48754035, -0.47483487,  0.20071602,  0.07649819]])

The variable 'ortho' is an orthogonal matrix as you can easily show by computing its determinant::

    In [6]: np.linalg.det(ortho)
    Out[6]: 1.0

Now, for a 7-dimension space we should expect 21 generalized Euler angles. We can get them by calling the function::

    In [8]: angles_list = pychemia.utils.mathematics.gea_all_angles(ortho)

    In [9]: np.array(angles_list)
    Out[9]:

    array([ -0.1392 , -0.03975, -0.37052, -0.48339,  0.13503,  0.39547, -0.22866, -0.31214, -0.57967,  0.84238,
            -2.66098,  0.05294,  0.3535 , -0.83619, -0.06953, -0.51353, -1.01991,  2.28529, -0.25373, -0.33108,
            -2.27736])

In fact we got 21 angles, all the angles in the range :math:`[-\pi, \pi]`
We can rebuild the original orthogonal matrix from those angles with::

    In [10]: matrix=pychemia.utils.mathematics.gea_orthogonal_from_angles(angles)

In fact, we can verify that we recover the original matrix::

    In [11]: np.max(matrix-ortho)
    Out[11]: 3.3306690738754696e-16

