#!/usr/bin/python -u

import os
import cPickle as pick


def main():
    if os.path.isfile('sdata.bin'):
        f = open('sdata.bin')
        evodata = pick.load(f)
        f.close()
    else:
        print 'Error: sdata'
        return 1

    x = evodata[1][1].saddle.get_e()
    print x


if __name__ == '__main__':
    main()
