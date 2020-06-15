worms = {'WormMeasGiw': 2, 
         'WormMeasGtau': 2, 
         'WormMeasGSigmaiw': 3,
         'WormMeasG4iw': 4,
         'WormMeasH4iw': 5,
         'WormMeasP2iwPH': 6,
         'WormMeasP2tauPH': 6,
         'WormMeasP2iwPP': 7,
         'WormMeasP2tauPP': 7,
         'WormMeasP3iwPH': 8,
         'WormMeasP3iwPP': 9,
         'WormMeasQQ': 10,
         'WormMeasQQtau': 10,
         'WormMeasQQQQ': 11,
         'WormMeasNQQdag': 12,
         'WormMeasQQdd': 13,
         'WormMeasUcaca': 14,
         'WormMeasUcacatau': 14,
         'WormMeasUccaa': 15,
         'WormMeasUccaatau': 15,
         'WormMeasQUDdag': 16}

def get_sector_index(cfg):
    """
    Get the worm sector number corresponding to the parameters WormMeas...
    """
    sectors = []
    for worm in worms.keys():
        if cfg[worm] != 0:
            sectors.append(worms[worm])
    sectors = list(set(sectors))
    if len(sectors) > 1:
        raise ValueError('Multiple worm sectors {} given.'.format(sectors))
    elif len(sectors) < 1:
        raise ValueError('No worm sector given')
    else:
        sector = sectors[0]

    if sector in [4, 5, 11] and cfg['FourPnt'] != 8:
        raise ValueError('For measuring three-frequency objects in worm sampling, set FourPnt = 8.')

    return sector
