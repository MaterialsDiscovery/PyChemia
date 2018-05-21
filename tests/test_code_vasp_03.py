import sys
import numpy
import pychemia

if pychemia.HAS_MAYAVI:

    from mayavi.mlab import quiver

    if not pychemia.HAS_MAYAVI:
        sys.exit(1)

    def test_quiver3d():
        x, y, z = numpy.mgrid[-2:3, -2:3, -2:3]
        r = numpy.sqrt(x ** 2 + y ** 2 + z ** 4)
        u = y * numpy.sin(r) / (r + 0.001)
        v = -x * numpy.sin(r) / (r + 0.001)
        w = numpy.zeros_like(z)
        obj = quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)
        return obj
