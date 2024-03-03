import numpy as np

class Trajectory:
    def __init__(
        self,
        f: np.array,
        g: np.array,
        terminal_time: float,
        initial_value: float = 1,
        dt: float = 10**-3,
    ):
        """
        Simulation of a diffusion process with Euler-Mauyama method,
                    dX(t) = f(t,X(t))dt + G(t,X(t))dW(t).
            - f: drift.
            - g: noise intensity,
            - terminal_time
            - initial_value: initial value. Set as 1 if none.
            - df: time step. Set as 10**-3 if none.

        """
        self.f = f
        self.g = g
        self.terminal_time = terminal_time
        self.initial_value = initial_value
        self.dt = dt

    def simulate(self):
        print("Start computing the path...")

        self.x = np.arange(start=0, stop=self.terminal_time, step=self.dt)
        self.y = np.zeros((np.shape(self.x)))
        self.y[0] = self.initial_value

        N = int(self.terminal_time/self.dt) + 1
        brownian_motion = np.cumsum(np.random.normal(loc = 0, scale = np.sqrt(self.dt), size = N))
        for k in range(1, len(self.x)):
            self.y[k] = (
                self.y[k - 1]
                + self.f * self.y[k - 1] * self.dt
                + self.y[k - 1] * self.dt * brownian_motion[k]
            )

        print("Computing the path is done...")
