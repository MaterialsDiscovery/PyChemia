import numpy as np


class Verlet:
    def __init__(self, x0, v0, delta, accel_function, dissipation=1.0, niter=None, min_vel=None):

        self.x0 = np.array(x0).reshape((-1, 3))
        self.v0 = np.array(v0).reshape((-1, 3))
        self.delta = delta
        self.accel_function = accel_function
        self.trajectory = [self.x0]
        self.niter = niter
        # Minimal velocity to stop
        self.min_vel = min_vel

    def run(self):

        x0 = self.x0
        x1 = self.x0 + self.v0 * self.delta + 0.5 * self.accel_function(self.x0) * self.delta ** 2
        self.trajectory.append(x1)
        iteration = 0
        while True:
            xn = 2 * x1 - x0 + self.accel_function(x1) * self.delta ** 2
            self.trajectory.append(xn)
            x0 = x1
            x1 = xn
            iteration += 1

            if self.niter is not None:
                if iteration == self.niter:
                    break

            if self.min_vel is not None:
                vel = np.array((x1 - x0) / self.delta)
                if np.max(np.apply_along_axis(np.linalg.norm, 1, vel)) > self.min_vel:
                    break
