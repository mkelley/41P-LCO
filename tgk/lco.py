# Licensed under a MIT style license - see LICENSE
"""lco - Las Cumbres Observatory info and parameters."""

# https://lco.global/observatory/mpccodes/
mpc_codes = {
    ('ogg', 'clma', '2m0a'): 'F65',
    ('ogg', 'clma', '0m4b'): 'T04',
    ('elp', 'doma', '1m0a'): 'V37',
    ('lsc', 'doma', '1m0a'): 'W85',
    ('lsc', 'domb', '1m0a'): 'W86',
    ('lsc', 'domc', '1m0a'): 'W87',
    ('lsc', 'aqwa', '0m4a'): 'W89',
    ('cpt', 'doma', '1m0a'): 'K91',
    ('cpt', 'domb', '1m0a'): 'K92',
    ('cpt', 'domc', '1m0a'): 'K93',
    ('coj', 'clma', '0m4b'): 'Q59',
    ('coj', 'doma', '1m0a'): 'Q63',
    ('coj', 'domb', '1m0a'): 'Q64',
    ('coj', 'clma', '2m0a'): 'E10',
    ('tfn', 'aqwa', '0m4a'): 'Z21',
    ('tlv', 'doma', '1m0a'): '097',
}

# LCO filter to PS1 catalog filter
filter2PS1 = {
    'rp': 'r',
    'w': 'r',
    'gp': 'g'
}
