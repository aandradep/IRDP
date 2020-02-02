import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import math
import datetime

class Bond:
    def __init__(self,
                 cupon: float = 10,
                 rate: float = 0.0532,
                 periodicity: int = 2,
                 date_val: datetime.date = datetime.date(2002, 12, 31),
                 date_end: datetime.date = datetime.date(2012, 1, 31)
                 ):
        self._cupon = cupon
        self._rate = rate
        self._periodicity = periodicity
        self._date_val = date_val
        self._date_end = date_end
        self.price_rate = self.BondPriceYield()
        self.modified_duration = self.ModifiedDuration()
        self.macauly_duration = self.MacaulyDuration()
        self._cupon_dates = self.CuponDates()

    @property
    def cupon(self):
        return self._cupon

    @property
    def rate(self):
        return self._rate

    @property
    def periodicity(self):
        return self._periodicity

    @property
    def date_val(self):
        return self._date_val

    @property
    def date_end(self):
        return self._date_end

    def CuponDates(self) -> list:
        """ Calculate the cupon dates of a bond.
        # Arguments:
            date_val {datetime.date}: date of pricing
            date_end {datetime.date}: date of maturity of the bond
            periodicity {int}: number of payments per year
        """
        ceiling_days = math.ceil((self.date_end - self.date_val).days / 365 * self.periodicity)
        days_add = np.arange(365 / self.periodicity, ceiling_days * 365 / self.periodicity,
                             step=365 / self.periodicity)[::-1] * -1
        dates_cupon = [self.date_end + datetime.timedelta(days) for days in days_add]
        dates_cupon.append(self.date_end)
        return dates_cupon

    def BondPriceYield(self) -> float:
        """ Calculate the price of a bond.
        # Arguments:
            cupon {float}: the cupon the bond pays for a 100 notional
            rate {rate}: the interest rate observed in the market
            date_val {datetime.date}: date of pricing
        """
        list_cupon_dates = self.CuponDates()
        cupons = np.repeat(self.cupon, len(list_cupon_dates))
        cupons[::-1][0] = cupons[::-1][0] + 100
        time_steps = [(date_cupon - self.date_val).days / 365 for date_cupon in list_cupon_dates]
        fd = np.exp(-np.array(time_steps) * self.rate)
        price = sum(fd * cupons)
        return price

    def MacaulyDuration(self) -> float:
        """ Calculate the macauly duration of a bond.
        # Arguments:
            cupon {float}: the cupon the bond pays for a 100 notional
            rate {rate}: the interest rate observed in the market
            date_val {datetime.date}: date of pricing
        """
        list_cupon_dates = self.CuponDates()
        cupons = np.repeat(self.cupon, len(list_cupon_dates))
        cupons[::-1][0] = cupons[::-1][0] + 100
        time_steps = [(date_cupon - self.date_val).days / 365 for date_cupon in list_cupon_dates]
        fd = np.exp(-np.array(time_steps) * self.rate)
        mac_duration = sum((fd * cupons * time_steps) / self.BondPriceYield())
        return (mac_duration)

    def ModifiedDuration(self) -> float:
        """ Calculate the macauly duration of a bond.
        # Arguments:
            periodicity {int}:  number of payments per year
            rate {rate}: the interest rate observed in the market
        """
        mac_duration = self.MacaulyDuration()
        mod_duration = mac_duration / (1 + self._rate / self._periodicity)
        return mod_duration


class InterestRateCurve(Bond):
    def __init__(self,
                 vec_term: list = [],
                 vec_rates: list = [],
                 vec_cupons: list = []
                 ):
        self._vec_term = vec_term
        self._vec_rates = vec_rates
        self._vec_cupons = vec_cupons

    def RawInterpolation(self, x_term: list = [], y_rate: list = [],
                         term_interpolate: float = 5) -> float:
        """ Interpolate the curve to get the rate for xout.
        # Arguments:
            x_term {list}:  vector with term for each rate in vec_rates
            y_rate {list}: vector with rates to use in the interpolation
            term_interpolate {float}: term to interpolate the curve and get an interest rate
        """
        next_position = np.searchsorted(x_term, term_interpolate)
        last_position = next_position - 1
        last_position = np.where(last_position == -1, None, last_position).tolist()
        next_position = np.where(next_position == len(x_term), None, next_position).tolist()

        if last_position is None:
            tao_1 = None
            tao_2 = x_term[next_position]
            rate_2 = y_rate[next_position]
        elif next_position is None:
            tao_1 = x_term[last_position]
            tao_2 = None
            rate_1 = y_rate[last_position]
        else:
            tao_1 = x_term[last_position]
            tao_2 = x_term[next_position]
            rate_1 = y_rate[last_position]
            rate_2 = y_rate[next_position]

        if (last_position is None) & (type(tao_2) == float):
            interpolation = rate_2
        elif (type(tao_1) == float) & (next_position is None):
            interpolation = rate_1
        elif (type(tao_1) == float) & (type(tao_2) == float):
            interpolation = ((term_interpolate - tao_1) / (tao_2 - tao_1) * (tao_2 / term_interpolate) * rate_2 +
                             (tao_2 - term_interpolate) / (tao_2 - tao_1) * (tao_1 / term_interpolate) * rate_1)
        return interpolation

    def _FunctionRate(self, Bond, x_term, y_rate) -> float:
        """ Get the rate of a bond based on a interest rate curve and its price.
        # Arguments:
            Bond {object}: an object of class Bond that has:
                - CuponDates() function to generate the dates of cupon
                - date_val {datetime.date} variables with the date of pricing
                - cupon {float} variable with the cupon of the bond
        # Detailes:
            The function leverage all the data in the Bond object for the operations. 
        """
        list_cupon_dates = Bond.CuponDates()
        time_steps = [(date_cupon - Bond.date_val).days / 365 for date_cupon in list_cupon_dates]
        rates_use = [self.RawInterpolation(x_term=x_term, y_rate=y_rate,
                                           term_interpolate=term) for term in time_steps]
        cupons = np.repeat(Bond.cupon, len(list_cupon_dates))
        N = len(list_cupon_dates)

        if N > 1:
            cupons_part = np.dot(np.array(cupons[0:(N - 1)]),
                                 np.exp(-np.array(rates_use[0:(N - 1)]) * np.array(time_steps[0:(N - 1)])))
        else:
            cupons_part = 0

        price_part = Bond.price_rate
        rate = (1 / time_steps[N - 1]) * (np.log(cupons[N - 1] + 100) - np.log(price_part - cupons_part))
        return rate

    def PriceCurve(self, Bond, x_term: list = [], y_rates: list = []) -> float:
        """ Price a bond based ona interest rate curve.
        # Arguments:
            Bond {object}: an object of class Bond that has:
                - CuponDates() function to generate the dates of cupon
                - date_val {datetime.date} variables with the date of pricing
                - cupon {float} variable with the cupon of the bond
        # Detailes:
            The function leverage all the data in the Bond object for the pricing. 
        """
        list_cupon_dates = Bond.CuponDates()
        time_steps = [(date_cupon - Bond.date_val).days / 365 for date_cupon in list_cupon_dates]
        rates_use = [self.RawInterpolation(x_term=x_term, y_rate=y_rates,
                                           term_interpolate=term) for term in time_steps]
        cupons = np.repeat(Bond.cupon, len(list_cupon_dates))
        cupons[::-1][0] = cupons[::-1][0] + 100
        price = np.dot(np.exp(-np.array(rates_use) * np.array(time_steps)), np.array(cupons))
        return price

    def CalibrateCurve(self, maturities: list = [], cupons: list = [],
                       date_calib: datetime.date = datetime.date.today()) -> np.array:
        bonds = [Bond(cupon=cupon, rate=rate, periodicity=1,
                      date_val=fechas_analisis, date_end=maturity) for
                 cupon, rate, maturity in zip(cupons, self._vec_rates, maturities)]
        dirty_price = [bond.price_rate for bond in bonds]

        x_term = self._vec_term.copy()
        y_term = self._vec_rates.copy()

        j = 0
        while j < 30:
            for i in np.arange(0, len(self._vec_rates)):
                y_term[i] = self._FunctionRate(bonds[i], x_term, y_term)
            j = j + 1

        curve = np.c_[x_term, y_term]
        price_curve = [self.PriceCurve(bond, x_term, y_term) for bond in bonds]
        dif_precios = sum(abs(np.array(dirty_price) - np.array(price_curve))) / len(price_curve)
        print(f'The absolute price difference is {round(dif_precios, 4)}')

        return curve

    def ForwardRate(self, t2: float = 3, t1: float = 2) -> float:
        """ Price a bond based ona interest rate curve.
        # Arguments:
            Bond {object}: an object of class Bond that has:
                - CuponDates() function to generate the dates of cupon
                - date_val {datetime.date} variables with the date of pricing
                - cupon {float} variable with the cupon of the bond
        # Detailes:
            The function leverage all the data in the Bond object for the pricing. 
        """
        denominator = t2 - t1
        numerator = (self.RawInterpolation(x_term=self._vec_term, y_rate=self._vec_rates, term_interpolate=t2) * t2 -
                     self.RawInterpolation(x_term=self._vec_term, y_rate=self._vec_rates, term_interpolate=t1) * t1)
        fwd_rate = numerator / denominator
        return fwd_rate


Y = [0.04531, 0.04839, 0.05425, 0.05990, 0.06440, 0.06684, 0.06920, 0.07055]
X = [0.5205479, 1.3890411, 3.1671233, 5.3917808, 7.4821918, 9.1561644, 11.5479452, 13.3315068]
maturities = [datetime.date(2019, 9, 11), datetime.date(2020, 7, 24),
              datetime.date(2022, 5, 4), datetime.date(2024, 7, 24), datetime.date(2026, 8, 26),
              datetime.date(2028, 4, 28), datetime.date(2030, 9, 18), datetime.date(2032, 6, 30)]
cupons = [7, 11, 7, 10, 7.5, 6, 7.75, 7]

fechas_analisis = datetime.date(2019, 3, 5)

initial_curve = InterestRateCurve(vec_term=X, vec_rates=Y)
curva = initial_curve.CalibrateCurve(maturities=maturities, cupons=cupons, date_calib=fechas_analisis)
curva[:, 1].tolist()
fig, ax = plt.subplots()
ax.scatter(curva[:, 0], curva[:, 1])
line = mlines.Line2D(X, Y, color='red')
ax.add_line(line)

test_bond = Bond(cupon=cupons[0], rate=Y[0], periodicity=1, date_val=fechas_analisis,
                 date_end=maturities[0])
test_bond.price_rate
test_bond.macauly_duration
test_bond.modified_duration

initial_curve.PriceCurve(test_bond, curva[:, 0].tolist(), curva[:, 1].tolist())
InterestRateCurve(vec_term=curva[:, 0].tolist(), vec_rates=curva[:, 1].tolist()).ForwardRate(t2=10, t1=4)
