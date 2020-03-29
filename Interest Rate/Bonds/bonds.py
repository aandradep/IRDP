import numpy as np
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


class InterpolationMethods:
    def __init__(self,
                 x_term: list = [],
                 y_rate: list = [],
                 term_interpolate: float = 5
            ):
        
        self._x_term = x_term
        self._y_rate = y_rate
        self._term_interpolate = term_interpolate    
    
    def RawInterpolation(self) -> float:
        """ Interpolate the curve to get the rate for xout.
        # Arguments:
            x_term {list}:  vector with term for each rate in vec_rates
            y_rate {list}: vector with rates to use in the interpolation
            term_interpolate {float}: term to interpolate the curve and get an interest rate
        """
        next_position = np.searchsorted(self._x_term, self._term_interpolate)
        last_position = next_position - 1
        last_position = np.where(last_position == -1, None, last_position).tolist()
        next_position = np.where(next_position == len(self._x_term), None, next_position).tolist()

        if last_position is None:
            tao_1 = None
            tao_2 = float(self._x_term[next_position])
            rate_2 = self._y_rate[next_position]
        elif next_position is None:
            tao_1 = float(self._x_term[last_position])
            tao_2 = None
            rate_1 = self._y_rate[last_position]
        else:
            tao_1 = float(self._x_term[last_position])
            tao_2 = float(self._x_term[next_position])
            rate_1 = self._y_rate[last_position]
            rate_2 = self._y_rate[next_position]

        if (last_position is None) & (type(tao_2) == float):
            interpolation = rate_2
        elif (type(tao_1) == float) & (next_position is None):
            interpolation = rate_1
        elif (type(tao_1) == float) & (type(tao_2) == float):
            interpolation = ((self._term_interpolate - tao_1) / (tao_2 - tao_1) * (tao_2 / self._term_interpolate) * rate_2 +
                             (tao_2 - self._term_interpolate) / (tao_2 - tao_1) * (tao_1 / self._term_interpolate) * rate_1)
        return interpolation

class InterestRateCurve(Bond,InterpolationMethods):
    def __init__(self,
                 vec_term: list = [],
                 vec_rates: list = []
                 ):
        self._vec_term = vec_term
        self._vec_rates = vec_rates
        
        
    def DiscountFactors(self, Bond) -> list:
        list_cupon_dates = Bond.CuponDates()
        time_steps = [(date_cupon - Bond.date_val).days / 365 for date_cupon in list_cupon_dates]
        rates_use = [InterpolationMethods(self._vec_term, self._vec_rates, term).RawInterpolation() 
                     for term in time_steps]
        df = np.exp(-np.array(rates_use) * np.array(time_steps))
        return df
        

    def PriceCurve(self, Bond) -> float:
        """ Price a bond based on a interest rate curve.
        # Arguments:
            Bond {object}: an object of class Bond that has:
                - CuponDates() function to generate the dates of cupon
                - date_val {datetime.date} variables with the date of pricing
                - cupon {float} variable with the cupon of the bond
        # Detailes:
            The function leverage all the data in the Bond object for the pricing. 
        """
        list_cupon_dates = Bond.CuponDates()
        df = self.DiscountFactors(Bond)
        cupons = np.repeat(Bond.cupon, len(list_cupon_dates))
        cupons[::-1][0] = cupons[::-1][0] + 100
        price = np.dot(df, np.array(cupons))
        return price

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
        numerator = (InterpolationMethods(self._vec_term,self._vec_rates,t2).RawInterpolation() * t2 -
                     InterpolationMethods(self._vec_term,self._vec_rates,t1).RawInterpolation() * t1)
        fwd_rate = numerator / denominator
        return fwd_rate
    
    
class CalibrateCurve(InterestRateCurve, Bond, InterpolationMethods):
    def __init__(self,
                 maturities: list = [],
                 vec_rates: list = [],
                 vec_term: list = [],
                 cupons: list = [],
                 date_calib: datetime.date = datetime.date.today()
                 ):
        self._maturities= maturities
        self._vec_rates = vec_rates
        self._vec_term = vec_term
        self._cupons = cupons
        self.date_calib = date_calib
    
    def _FunctionRate(self, Bond) -> float:
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
        rates_use = [InterpolationMethods(self._vec_term , self._vec_rates ,term).RawInterpolation()
                     for term in time_steps]
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

    def Calibration(self) -> np.array:
        bonds = [Bond(cupon=cupon, rate=rate, periodicity=1,
                      date_val=self.date_calib, date_end=maturity) for
                 cupon, rate, maturity in zip(self._cupons, self._vec_rates, self._maturities)]
        dirty_price = [bond.price_rate for bond in bonds]
        
        x_term = self._vec_term.copy()
        y_term = self._vec_rates.copy()
        
        j = 0
        while j < 30:
            for i in np.arange(0, len(self._vec_rates)):
                y_term[i] = self._FunctionRate(bonds[i])
            j = j + 1
    
        curve = np.c_[x_term, y_term]
        inter_curva = InterestRateCurve(x_term,y_term)
        price_curve = [inter_curva.PriceCurve(bond) for bond in bonds]
        dif_precios = sum(abs(np.array(dirty_price) - np.array(price_curve))) / len(price_curve)
        print(f'The absolute price difference is {round(dif_precios, 4)}')
        return curve
    