import astropy
from astropy.table import Table, Column
import numpy as np
import math


class Formido:
    def __init__(self):
        self.corrections = {"distmod": 0, "ebv": 0, "ebv_r": 0, "ebv_b": 0}
        self.stars = Table()

    # Generates the needed magnitude adjustments to adjust generated isochrones.
    def Adjustments(self, dist=0, ebv=0, rb_scalar=0, bb_scalar=0):
        self.corrections["distmod"] = 5.0 * np.log10(dist) - 5
        self.corrections["ebv"] = ebv
        self.corrections["ebv_r"] = ebv * rb_scalar
        self.corrections["ebv_b"] = ebv * bb_scalar

        return {
          'd_rmag': self.corrections["distmod"] + (ebv * rb_scalar),
          'd_bmag': self.corrections["distmod"] + (ebv * bb_scalar),
          'd_color': (self.corrections["distmod"] + (ebv * rb_scalar)) - (self.corrections["distmod"] + (ebv * bb_scalar))
        }

    # Adds stars and their associated bands to construct CMDs.
    def AddStars(self, ra, dec, red_band=[], blue_band=[]):
        if len(ra) == len(dec):
            self.stars["ra"] = ra
            self.stars["dec"] = dec
        else:
            raise Exception(
                "Coordinates provided are inconsistent, arrays of different sizes provided."
            )

        if red_band == []:
            self.stars["rb"] = np.zeros(len(ra))
        else:
            self.stars["rb"] = red_band

        if blue_band == []:
            self.stars["bb"] = np.zeros(len(ra))
        else:
            self.stars["bb"] = blue_band

    # Subtracts the bluer band by the redder band.
    def GenerateColor(self):
        blue = self.stars["bb"]
        red = self.stars["rb"]
        return Column(data=(blue - red), name="b_min_r")

    # Generates a probability based on the distance relative to the half-light radius.
    def GenerateHLProb(self, eff_rad, center={"ra": 0, "dec": 0}, eps=0, theta=0):
        ra = self.stars["ra"]
        dec = self.stars["dec"]
        prob = []

        for k in self.stars:
            r = np.sqrt(
                ((k["ra"] - center["ra"]) * np.cos(math.radians(k["dec"]))) ** 2
                + (k["dec"] - center["dec"]) ** 2
            )

            rh = eff_rad * (1 - eps) * (1 + eps * np.cos(theta))
            p_dist = np.exp(-(r ** 2) / (2 * (rh ** 2)))
            prob.append(p_dist)

        return Column(data=prob, name="hl_prob")

    # Generates probabilities based on the distance to an isochrone.
    def GenerateIDProb(self, iso_b, iso_r, sigma=0.1, verbose=False):
        iso_dist = []

        if len(iso_r) != len(iso_b):
            raise Exception("Inconsistent isochrone parameters, column size mismatch.")

        for star in np.arange(len(self.stars)):
            if star == len(self.stars):
                break
            else:
                color_diff = Column(
                    data=(self.stars[star]["bb"] - self.stars[star]["rb"])
                    - (iso_b - iso_r),
                    name="color_diff"
                )

                mag_diff = Column(data=self.stars[star]["rb"] - iso_r, name="mag_diff")
                
                dist_vector = np.sqrt((color_diff.data ** 2) + (mag_diff.data ** 2))
               
                min_dist = np.min(dist_vector)
                
                iso_dist.append(min_dist)
        iso_dist_col = Column(data=iso_dist, name="min_iso_dist")

        if verbose:

            print(iso_dist_col, np.min(iso_dist_col), np.max(iso_dist_col))

        idprob = Column(
            data=np.exp(-(iso_dist_col.data ** 2) / (2 * (sigma ** 2))), name="id_prob"
        )
        return idprob

        return Column(data=idprob, name="iso_prob")

        # star_dist = []
        # for ref in np.arange(len(iso_r)):
        #   if (ref == len(iso_r)):
        #     break
        #   else:
        #     blue = self.stars[star]['bb']
        #     red = self.stars[star]['rb']
        #     delta_br = (self.stars[star]['bb'] - self.stars[star]['rb']) - (iso_b[] - iso_r[ref])
        #     delta_r = self.stars[star]['rb'] - iso_r

    # Returns simple star information.
    def get(self):
        return self.stars

