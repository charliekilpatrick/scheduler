import Constants as C
from Target import TargetType, Target
import Logs

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np, csv, sys
from astropy.time import Time

class Telescope(metaclass=ABCMeta):

    @abstractmethod
    def set_targets(self, targets):
        pass

    @abstractmethod
    def get_targets(self):
        pass

    @abstractmethod
    def compute_exposures(self):
        pass

    @abstractmethod
    def write_schedule(self, observatory_name, obs_date, good_targets,
        output_files=None, fieldcenters=None, pointing=None):
        pass

    def round_to_num(self, round_to_num, input_to_round):
        return int(round_to_num*round(float(input_to_round)/round_to_num))

    def time_to_S_N(self, desired_s_n, apparent_mag, zeropoint, min_exp=45.):
        term1 = (desired_s_n / 3.)**2
        term2 = 0.4*(apparent_mag - zeropoint)
        exp_time = term1*10**term2

        if exp_time < min_exp:
            exp_time = min_exp

        return exp_time

    def write_fieldcenters(self, targets, fieldcenters_file, ampl='1'):

        with open(fieldcenters_file, 'w') as fc:

            fc.write('# field ampl ra dec epoch raD decD RAoffset DecOffset \n')

            for t in targets:

                name = t.name
                ra = t.coord.ra.degree
                dec = t.coord.dec.degree

                hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
                ra_hms = hmsdms.split()[0]
                dec_dms = hmsdms.split()[1]

                fc_line = f'{name:<40} {ampl} {ra_hms} {dec_dms} J2000  '+\
                          f'{ra:>11}  {dec:>11}    0.0000000    0.0000000 \n'

                fc.write(fc_line)

    def write_phot_file(self, targets, out_phot_file):

        with open(out_phot_file, 'w') as phot:

            for t in targets:

                name = t.name
                ra = t.coord.ra.degree
                dec = t.coord.dec.degree

                for filt in t.exposures.keys():

                    zeropoint = self.filters[filt]
                    exptime = t.exposures[filt]
                    m3sigma = self.limiting_magnitude(zeropoint, exptime, 3)
                    priority = t.orig_priority

                    phot_line = f'{name} {ra} {dec} {filt} {exptime} {m3sigma} {priority} \n'

                    phot.write(phot_line)

    def limiting_magnitude(self, zeropoint, exptime, s_n):
        return(2.5 * np.log10(exptime / (s_n/3.)**2) + zeropoint)

    def compute_net_priorities(self):

        targets = self.get_targets()

        total_p = np.sum([t.priority for t in targets])
        print("Total Priority: %s" % total_p)

        total_good_time = np.sum([t.total_observable_min for t in targets])
        print("Total Good Time: %s" % total_good_time)

        total_exp_time = np.sum([t.total_minutes for t in targets])
        print("Total Exposure Time: %s" % total_exp_time)

        total_prob = 0
        if (total_p > 0 and total_good_time > 0 and total_exp_time > 0):
            for t in targets:
                frac_p = float(t.priority) / float(total_p)
                frac_time = float(t.total_observable_min)/float(total_good_time)
                frac_exp_time = (1.0-float(t.total_minutes)/float(total_exp_time))

                if (frac_exp_time == 0.0):
                    frac_exp_time = 1.0

                total_prob += frac_p*frac_time*frac_exp_time

            for t in targets:

                frac_p = float(t.priority) / float(total_p)
                frac_time = float(t.total_observable_min)/float(total_good_time)
                frac_exp_time = (1.0-float(t.total_minutes)/float(total_exp_time))

                t.net_priority = t.priority+((frac_p*frac_time*frac_exp_time)/total_prob)
                print("Nat: %s; Net: %0.5f" % (t.priority, t.net_priority))
        else:
            print("No valid targets...")
