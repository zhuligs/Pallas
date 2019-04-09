import itin


def opt(amode):
    if itin.interface == 'lammps':
        reac = oljobsw()
    return reac