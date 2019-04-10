'''
Find parameters with SpaceEx
'''

import time

# make sure hybridpy is on your PYTHONPATH: hyst/src/hybridpy
import hybridpy.hypy as hypy

def main():
    '''main entry point for vanderpol_althoff script'''

    bottom = 0
    top = 2

    while True:
        mid = (top + bottom) / 2

        print "Trying [0, {}]".format(mid)

        init = 'abortmin==120 & abortmax==120 & t==0 & vx==0 & vy==0 & '.format(mid) + \
                      '-925<=x<=-875 & -425<=y<=-375 & Tmax==300 & loc()==P2"'

        #e = hypy.Engine('spaceex', '-aggregation none -output-format none')
        e = hypy.Engine('spaceex', '-output-format none')

        e.set_verbose(True)
        e.set_input('abort.xml')
        e.set_init(init)
        e.set_output('out.xml')

        result = e.run(parse_output=True)
        code = result['code']

        if code != hypy.Engine.SUCCESS:
            raise RuntimeError('Error: ' + str(code))

        if result['output']['safe'] != True:
            print "System was unsafe with abortmax:{}".format(mid)
            top = mid
        else:
            print "System was safe with abortmax:{}".format(mid)
            bottom = mid

        if top == bottom + 1:
            print "maximum safe range: [0, {}]".format(bottom)
            break
        
if __name__ == '__main__':
    main()









