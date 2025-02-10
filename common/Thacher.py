import common.Constants as C
from common.Target import TargetType, Target
import common.Logs
import common.Telescope
import common.Utilities

import numpy as np
import itertools

class Thacher(common.Telescope.Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Thacher"
        # Filter name: Zero-point
        self.filters = {
            C.V_band:14.93966783,
            C.r_band:15.1197765,
            C.i_band:14.52342636,
            C.z_band:13.68546838,
            C.g_band:15.29722848
        }
        self.exp_funcs = {
            TargetType.Supernova: self.compute_sn_exposure,
            TargetType.Template: self.compute_template_exposure,
            TargetType.Standard: self.compute_standard_exposure,
            TargetType.GW: self.compute_gw_exposure
        }

    def set_targets(self, targets):
        self.targets = targets

    def get_targets(self):
        return self.targets

    def compute_sn_exposure(self, sn):
        exposures = {}

        # Compute the current guess at apparent magnitude
        if sn.obs_date is None or sn.ref_date is None:
            days_from_disc = 0.
        else:
            days_from_disc = (sn.obs_date - sn.ref_date).jd

        mag_reduction = days_from_disc*0.03
        adj_app_mag = sn.apparent_mag + mag_reduction

        # Change S/N depending on phase...
        s_to_n = 10 # base signal to noise
        if days_from_disc <= 10:
            s_to_n = 30
        elif days_from_disc > 10 and days_from_disc <= 60:
            s_to_n = 20

        exposures.update({C.r_band: self.round_to_num(C.round_to, 
            self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.r_band]))})
        exposures.update({C.i_band: self.round_to_num(C.round_to, 
            self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.i_band]))})

        # Finally, don't go less than 45s (~ readout time), don't go more than 600s on Swope
        for key, value in exposures.items():
            if exposures[key] < 45:
                exposures[key] = 45
            elif exposures[key] > 600:
                exposures[key] = 600

        sn.exposures = exposures

    def compute_standard_exposure(self, std):
        exposures = {}
        s_to_n = 100

        # Don't know what the apparent mag should be?
        exposures.update({C.r_band: self.round_to_num(C.round_to, 
            self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.r_band]))})
        exposures.update({C.i_band: self.round_to_num(C.round_to, 
            self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.i_band]))})
        exposures.update({C.B_band: self.round_to_num(C.round_to, 
            self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.B_band]))})
        exposures.update({C.V_band: self.round_to_num(C.round_to, 
            self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.V_band]))})

        # Finally, don't go less than 10s for Nickel std, don't go more than 600s on Nickel
        for key, value in exposures.items():

            if exposures[key] < 10:
                exposures[key] = 10
            elif exposures[key] > 600:
                exposures[key] = 600

        std.exposures = exposures

    def compute_template_exposure(self, tmp):
        exposures = {}
        exposures.update({C.B_band: 1800})
        exposures.update({C.V_band: 1200})
        exposures.update({C.r_band: 1200})
        exposures.update({C.i_band: 1200})

        tmp.exposures = exposures


    def compute_gw_exposure(self, gw):
        exposures = {}

        s_to_n = 5.
        i_exp = self.time_to_S_N(s_to_n, gw.apparent_mag, self.filters[C.i_band])
        mean_exp = self.round_to_num(C.round_to, i_exp)

        exposures.update({C.i_band: mean_exp})

        gw.exposures = exposures

    def compute_exposures(self):
        for tgt in self.targets:

            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])

            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to/60.

                self.exp_funcs[tgt.type](tgt) # Sets exposures for each target by target type

                # per observatory
                fudge_factor = 120. if len(tgt.exposures) > 0 else 0 # Build in a fudge factor based on # of exps GW

                total_minutes = (sum(tgt.exposures.values()) + fudge_factor)/60. # Sum total minutes
                tgt.total_minutes = round(total_minutes * 60./C.round_to) * C.round_to/60.
                only_exposures = sum(tgt.exposures.values())/60.
                tgt.total_minutes_only_exposures = round(only_exposures * 60./C.round_to) * C.round_to/60.

            good_airmass = tgt.raw_airmass_array[np.asarray(np.where(tgt.raw_airmass_array <= C.airmass_threshold))]
            integrated_good_am = np.sum(good_airmass)
            if integrated_good_am > 0:
                tgt.total_good_air_mass = integrated_good_am

            if tgt.total_observable_min > 0:
                tgt.fraction_time_obs = float(tgt.total_minutes)/float(tgt.total_observable_min)

    def filter_row(self, exp_name, exp_time):
        filter_row = []
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(exp_name)
        filter_row.append(exp_time)

        return filter_row

    def make_observation_dict(self, name,
        observer = 'Charlie Kilpatrick', max_airmass = 2.5, priority = 50,
        num_images = '1', interval = '0', iterate = '1', focus_on = 'False',
        bestefforts = 'True', pa = 0, binvalue = 1,
        project_name = 'Thacher Gravitational Wave Astronomy (TGWA)',
        description = 'Gravitational Wave follow up by Thacher Observatory'):

        default = {
            'Request': {
                'bestefforts': bestefforts,
                'ID': {
                    'name': name
                },
                'UserName': {
                    'name': observer
                },
                'Observers': {
                    'name': observer
                },
                'Description': {
                    'name': description
                },
                'Reason': {
                    'name': 'Monitor=1'
                },
                'Project': {
                    'name': project_name
                },
                'Schedule': {
                    'Airmass': {
                        'name': max_airmass
                    },
                    'Priority': {
                        'name': priority
                    }
                }
            }
        }

        return(default)

    def make_target_dict(self, objname, coord, filt, exptime,
        observer = 'Charlie Kilpatrick', max_airmass = 2.5, priority = 1000,
        num_images = '1', interval = '0', iterate = '1', focus_on = 'False',
        bestefforts = 'True', pa = 0, binvalue = 1,
        project_name = 'Thacher Gravitational Wave Astronomy (TGWA)',
        description = 'Gravitational Wave follow up by Thacher Observatory'):

        ra_deg = coord.ra.degree
        dec_deg = coord.dec.degree

        default = {
            'Target': {
                    'count': num_images,
                    'interval': interval,
                    'autofocus': focus_on,
                    'Name': {
                        'name': objname
                    },
                    'Description': {
                        'name': description
                    },
                    'Coordinates': {
                        'RightAscension': {
                            'name': ra_deg
                        },
                        'Declination': {
                            'name': dec_deg
                        }
                    },
                    'PositionAngle': {
                        'name': pa
                    },
                    'Picture': {
                        'count': iterate,
                        'Name': {
                            'name': filt
                        },
                        'ExposureTime': {
                            'name': exptime
                        },
                        'Binning': {
                            'name': binvalue
                        },
                        'Filter': {
                            'name': filt
                        }
                    }
                }
            }

        return(default)

    def make_header_dict(self, contact = 'Charlie Kilpatrick',
        email = 'ckilpatrick@northwestern.edu', xml_version = '1.0',
        encoding = 'ISO-8859-1', RTML_version = '2.3'):
        default = {
            'RTML': {
                'version': RTML_version,
                'Contact': {
                    'User': {
                        'name': contact
                    },
                    'Email': {
                        'name': email
                    }
                }
            }
        }

        return(default)

    def serialize_xml(self, targets, outfile):

        xml_header = self.make_header_dict()

        f = open(outfile, 'w')
        i = 0

        # Make an observation Request object
        for tgt in targets:
            name = tgt.name
            coord = tgt.coord

            for filt in tgt.exposures.keys():

                number = str(i).zfill(7)
                exptime = str(tgt.exposures[filt])  
                filt = str(filt)

                obs = self.make_observation_dict(name, priority = 100,
                    project_name='Thacher Gravitational Wave Observing')

                xml_header['RTML']['Request'+number] = obs['Request']

                target = self.make_target_dict(tgt.name.lower(), coord, filt, 
                    str(tgt.exposures[filt]))

                xml_header['RTML']['Request'+number]['Target']=target['Target']

                i += 1

        # This is the order that the keys need to have in the final
        # xml file.  Dictionaries don't preserve order whereas the
        # Thacher software cares about the order of each key.
        # The dicttoxml function will write out keys for each layer
        # in the order specified here.
        order = ['xml', 'version', 'encoding', 'RTML', 'Contact', 'User',
                     'Email', 'Request', 'ID' , 'UserName',
                     'Observers', 'Name', 'Description', 'Reason',
                     'Project', 'Schedule', 'Airmass', 'Priority',
                     'Target', 'count', 'interval', 'autofocus',
                     'Coordinates', 'RightAscension', 'Declination',
                     'Picture', 'ExposureTime', 'Binning', 'Filter']
        common.Utilities.dicttoxml(f, xml_header, depth = 0, name_name='name', 
            order=order)

        f.close()

    def write_schedule(self, observatory_name, obs_date, targets, outdir=None,
        output_files=None, fieldcenters=None, pointing=None, **kwargs):

        if not output_files:
            output_files = "%s_%s_Schedule" % (self.name, obs_date.strftime('%Y%m%d'))

        output_rows = self.serialize_output_rows(targets, pointing=pointing)

        if fieldcenters:
            self.write_fieldcenters_file(targets, fieldcenters)

        self.write_csv_output(output_rows, output_files+'.csv')
        self.write_phot_file(targets, output_files+'.phot')
        self.serialize_xml(targets, output_files+'.xml')
