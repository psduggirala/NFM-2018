'''
emsoft-2019 SpaceEx measurements (generates plot/*.dat files for plotting)
'''

import time

# make sure hybridpy is on your PYTHONPATH: hyst/src/hybridpy
import hybridpy.hypy as hypy

def spaceex(abortmin, abortmax, agg="none", tol=0.1, directions='box', timeout=5):
    '''
    run spaceex with the given settings. will run multiple times if time is small
    '''

    small = 10.0

    rv = spaceex_single(abortmin, abortmax, agg=agg, tol=tol, directions=directions, timeout=timeout)

    if rv[0] == 'safe' and rv[1] < small:
        # run four more times and return average of middle three runs
        times = [rv[1]]

        for _ in range(4):
            rv = spaceex_single(abortmin, abortmax, agg=agg, tol=tol, directions=directions, timeout=timeout)
            times.append(rv[1])

        times.sort()

        rv[1] = sum(times[1:4]) / 3.0

        print "Runtimes: {}, returning {}".format(times, rv[1])

    return rv

def spaceex_single(abortmin, abortmax, agg="none", tol=0.1, directions='box', timeout=5):
    '''
    run spaceex with the given settings. runs only one time

    Switching to passive occurs at time [abortmin, abortmax]

    returns ('safe', runtime_sec), ('unsafe', runtime_sec)
    '''

    init = 'abortmin=={} & abortmax=={} & t==0 & vx==0 & vy==0 & '.format(abortmin, abortmax) + \
              '-925<=x<=-875 & -425<=y<=-375 & Tmax==250 & loc()==P2"'

    params = '-output-format none -aggregation {} -flowpipe_tol {} -directions {}'.format(
        agg, tol, directions)
    e = hypy.Engine('spaceex', params)

    e.set_verbose(True)
    e.set_input('abort.xml')
    e.set_init(init)
    e.set_output('out.xml')

    print "Running Spaceex with r={} and params {}...".format(70-abortmin, params),
    result = e.run(parse_output=True, print_stdout=False, timeout=2*timeout) # use 2*timeout as a hard timeout
    code = result['code']

    if code == hypy.Engine.TIMEOUT_TOOL:
        rv = ['unsafe', None] # no time means it was a timeout
    elif code != hypy.Engine.SUCCESS:
        raise RuntimeError('Error: ' + str(code))
    elif result['output']['safe']:
        rv = ['safe', result['tool_time']]
    else:
        rv = ['unsafe', result['tool_time']]

    print rv

    return rv

def linear_search(f, safe_rad, unsafe_rad, agg='chull', tol=0.1, directions='box', timeout=5):
    'do a linear seach up to the unsafe value, printing the results to f'

    for rad in range(safe_rad + 1, unsafe_rad):
        result = spaceex(70-rad, 70+rad, agg=agg, tol=tol, directions=directions, timeout=timeout)

        result, sec = result # extract result tuple

        if result == 'unsafe':
            break

        f.write("{} {}\n".format(rad, sec))

        if sec > timeout + 5:
            break

def run_chull(timeout):
    'run convex hull measurements, saves to plot/data_chull.dat'

    print "Running chull measurements"
    start = time.time()
    step = 5
    params = [(0.2, 'box'), (0.01, 'box'), (0.2, 'oct'), (0.01, 'oct')]


    with open('plot/data_chull.dat', 'w') as f:

        for param in params:
            tol, direction = param

            f.write('" {} (fp-tol={})"\n'.format(direction, tol))

            for rad in range(0, 75, step):
                result = spaceex(70-rad, 70+rad, agg='chull', tol=tol, directions=direction, timeout=timeout)

                result, sec = result # extract result tuple

                if result == 'unsafe':
                    # search with a smaller time step
                    linear_search(f, rad-5, rad, agg='chull', tol=tol, directions=direction, timeout=timeout)
                    break

                f.write("{} {}\n".format(rad, sec))

                if sec > timeout + 5: # don't need binary search here because we have a valid measurement
                    break

            f.write("\n\n")
            print "" # newline

    print "Finished chull in {} minutes".format((time.time() - start) / 60.0)

def run_unaggregated(timeout):
    'run unaggregated measurements, saves to plot/data_unagg.dat'

    print "Running unagg measurements"
    start = time.time()
    step = 5
    params = [(0.2, 'box'), (0.02, 'box'), (0.5, 'oct'), (0.2, 'oct')]

    with open('plot/data_unagg.dat', 'w') as f:

        for param in params:
            tol, direction = param

            f.write('" {} (fp-tol={})"\n'.format(direction, tol))

            for rad in range(0, 75, step):
                result = spaceex(70-rad, 70+rad, agg='none', tol=tol, directions=direction, timeout=timeout)

                result, sec = result # extract result tuple

                if result == 'unsafe':
                    # search with a smaller time step
                    linear_search(f, rad-5, rad, agg='none', tol=tol, directions=direction, timeout=timeout)
                    break

                f.write("{} {}\n".format(rad, sec))

                if sec > timeout: # don't need binary search here because we have a valid measurement
                    break

            f.write("\n\n")
            print "" # newline

    print "Finished unagg in {} minutes".format((time.time() - start) / 60.0)

def main():
    '''main entry point for vanderpol_althoff script'''

    # r= 56 is unsafe with tol = 0.1 & 0.05, safe with tol=0.025 (~30 sec),
    # also safe with tol=0.1 and dir='oct', but takes longer (~70 sec)
    # safe with tol=0.2 and dir='oct' (~42 sec)

    # chull/oct is safe up to 25

    timeout=60
    
    run_chull(timeout=timeout) # takes about 11 minutes with timeout=60
    
    run_unaggregated(timeout=timeout) # takes about 24 minutes with timeout=60
    return
        
if __name__ == '__main__':
    main()









