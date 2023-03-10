import configparser
import sys

config = configparser.ConfigParser()
config['MandatoryParams'] = {'ra': '0.0',
                             'dec': '0.0',
                             'ebv': '0.0271',
                             'err_ebv': '0.0005',
                             'aj': '0.014',
                             'ah': '0.009',
                             'ak': '0.006',
                             'PS1_DR': 'Dr2',
                             'Gaia_DR': 'Dr3'}
with open('config.ini', 'w') as configfile:
    config.write(configfile)
