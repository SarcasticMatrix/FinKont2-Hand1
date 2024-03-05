import numpy as np
from scipy.stats import norm


class BlackScholes:
    def __init__(
        self,
        T: float,
        strike: float,
        r: float,
        sigma: float,
        q: float = 0,
        option_type: str = "call",
        greek_type: str = "none",
    ):
        """
        Compute several formulas
        - T: maturity date.
        - strike: the strike.
        - r: interest rate.
        - sigma: volatlity.
        - q: dividends. Is set at 0 by default.
        - option_type: "call" or "put".
        - greek_type: none, delta, gamma

        - price: arbitrage free price
        - d1: d1
        - d2: d2
        """

        self.T = T
        self.strike = strike
        self.r = r
        self.sigma = sigma
        self.q = q

        if option_type in ["put", "call"]:
            self.option_type = option_type
        else:
            print("option_type must be: 'call' or 'put' ")

        if greek_type in ["none", "delta", "gamma"]:
            self.greek_type = greek_type
        else:
            print("greek_type must be: 'none', 'delta' or 'gamma' ")

    def call(self,spot):
        """ Compute call """
        cdf_value1 = norm().cdf(self.d1)
        cdf_value2 = norm().cdf(self.d2)
        self.price = (
            spot * np.exp(-self.q * self.T) * cdf_value1
            - self.strike * np.exp(-self.r * self.T) * cdf_value2
        )

    def put(self,spot):
        """ Compute put """
        cdf_value1 = norm().cdf(self.d1)
        cdf_value2 = norm().cdf(self.d2)
        self.price = (
            spot * np.exp(-self.q * self.T) * cdf_value1
            + self.strike * np.exp(-self.r * self.T) * cdf_value2
        )

    def call_delta(self):
        cdf_value1 = norm().cdf(self.d1)
        self.price = cdf_value1

    def put_delta(self):
        cdf_value1 = norm().cdf(self.d1)
        self.price = np.exp(-self.q * self.T) * (cdf_value1 - 1) # cdf_value1 + np.exp(-self.q * self.T)

    def valuation(self,spot):
        """
        Arbitrage free price of the option
        """

        self.d1 = (
            np.log(spot / self.strike)
            + ((self.r - self.q) + 0.5 * self.sigma**2) * self.T
        ) / (self.sigma * np.sqrt(self.T))
        self.d2 = self.d1 - self.sigma * np.sqrt(self.T)

        if self.greek_type == "none":
            if self.option_type == "call":
                self.call(spot)
            elif self.option_type == "put":
                self.put(spot)
        elif self.greek_type == "delta":
            if self.option_type == "call":
                self.call_delta()
            elif self.option_type == "put":
                self.put_delta()

        return self.price
