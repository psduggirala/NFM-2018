'''
Manually created hybrid automaton for hylaa testing.
'''

import math, time

from hylaa.timerutil import Timers

from hylaa.hybrid_automaton import LinearHybridAutomaton, LinearConstraint
from hylaa.hybrid_automaton import HyperRectangle
from hylaa.engine import HylaaSettings
from hylaa.engine import HylaaEngine
from hylaa.plotutil import PlotSettings

from numpy import array as nparray

def define_ha(r):
    '''make the hybrid automaton and return it'''

    # deagg: 120-140 = singular matrix ?? (invert used for aggregation)
    # deagg: 110-140 does not contain 140,140 = incorrect result ??
    # agg: [130, 140] = singular matrix
    # agg: [140, 140] doesn't do aggregation... only single set in successor (looks like 3 from the plot)

    #r = 65
    #abortmintime = 70 - r
    #abortmaxtime = 70 + r

    abortmintime = 70-r
    abortmaxtime = 70+r

    ha = LinearHybridAutomaton('Spacecraft Rendezvous')
    ha.variables = ["x", "y", "vx", "vy", "t", "affine"]

    p2 = ha.new_mode('P2')
    p2.a_matrix = nparray([\
        [0.0, 0.0, 1.0, 0.0, 0.0, 0], \
        [0.0, 0.0, 0.0, 1.0, 0.0, 0], \
        [-0.057599765881773, 0.000200959896519766, -2.89995083970656, 0.00877200894463775, 0.0, 0], \
        [-0.000174031357370456, -0.0665123984901026, -0.00875351105536225, -2.90300269286856, 0.0, 0.0], \
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], \
        [0, 0, 0, 0, 0, 0]])
    p2.c_vector = nparray([0.0, 0.0, 0.0, 0.0, 0.0, 0])

    inv1_1 = LinearConstraint([0.0, 0.0, 0.0, 0.0, 1.0, 0], abortmaxtime)
    inv1_2 = LinearConstraint([1.0, 0.0, 0.0, 0.0, 0.0, 0], -100.0)
    p2.inv_list = [inv1_1, inv1_2]

    p3 = ha.new_mode('P3')
    p3.a_matrix = nparray([\
        [0.0, 0.0, 1.0, 0.0, 0.0, 0], \
        [0.0, 0.0, 0.0, 1.0, 0.0, 0], \
        [-0.575999943070835, 0.000262486079431672, -19.2299795908647, 0.00876275931760007, 0.0, 0], \
        [-0.000262486080737868, -0.575999940191886, -0.00876276068239993, -19.2299765959399, 0.0, 0], \
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.],\
        [0, 0, 0, 0, 0, 0]])
    p3.c_vector = nparray([0.0, 0.0, 0.0, 0.0, 0.0, 0])

    inv2_1 = LinearConstraint([-1.0, 0., 0., 0., 0., 0], 100.0)
    inv2_2 = LinearConstraint([1.0, 0., 0., 0., 0., 0], 100.0)
    inv2_3 = LinearConstraint([0., -1.0, 0., 0., 0., 0], 100)
    inv2_4 = LinearConstraint([0., 1.0, 0., 0., 0., 0], 100)
    inv2_5 = LinearConstraint([-1.0, -1.0, 0.0, 0.0, 0.0, 0], 141.1)
    inv2_6 = LinearConstraint([1.0, 1.0, 0.0, 0.0, 0.0, 0], 141.1)
    inv2_7 = LinearConstraint([-1.0, 1.0, 0., 0., 0., 0], 141.1)
    inv2_8 = LinearConstraint([1.0, -1.0, 0., 0., 0., 0], 141.1)
    inv2_9 = LinearConstraint([0., 0., 0., 0., 1., 0], abortmaxtime)
    p3.inv_list = [inv2_1, inv2_2, inv2_3, inv2_4, inv2_5, inv2_6, inv2_7, inv2_8, inv2_9]

    passive = ha.new_mode('passive')
    passive.a_matrix = nparray(\
        [[0, 0, 1, 0, 0, 0], \
         [0, 0, 0, 1, 0, 0], \
         [0.0000575894721132000, 0, 0, 0.00876276, 0, 0], \
         [0, 0, -0.00876276, 0, 0, 0], \
         [0, 0, 0, 0, 0, 1.], \
         [0, 0, 0, 0, 0, 0]])

    # x'==vx & y'==vy & vx'==0.0000575894721132000*x+0.00876276*vy & vy'==-0.00876276*vx & t'==1

    passive.c_vector = nparray([0, 0, 0, 0, 0., 0])
    #inv3_1 = LinearConstraint([0.0, 0.0, 0.0, 0.0, -1.0], -120)
    #passive.inv_list = [inv3_1]

    error = ha.new_mode('error')
    error.is_error = True

    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, -1., 0., 0., 0., 0], 0)]

    # In the aborting mode, the vehicle must avoid the target, which is modeled as a box B with
    # 0.2m edge length and the center placed as the origin
    rad = 5
    ha.new_transition(passive, error).condition_list = [ \
        LinearConstraint([1, 0, 0., 0., 0., 0], rad), \
        LinearConstraint([-1, 0, 0., 0., 0., 0], rad), \
        LinearConstraint([0, 1., 0., 0., 0., 0], rad), \
        LinearConstraint([0, -1., 0., 0., 0., 0], rad)]

    #In the rendezvous attempt the spacecraft must remain within the lineof-sight
    #cone L = {[x, y]^T | (x >= -100m) AND (y >= x*tan(30)) AND (-y >= x*tan(30))}
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([1, 0, 0., 0., 0., 0], -100)]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([-0.57735, 1, 0, 0., 0., 0], 0)]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([-0.57735, -1, 0., 0., 0., 0], 0)]

    # sqrt(vx^2 + vy^2) should stay below 0.055 m/SECOND (time in model is in MINUTES)
    meters_per_sec_limit = 0.055
    meters_per_min_limit = meters_per_sec_limit * 60
    h = meters_per_min_limit * math.cos(math.pi / 8.0)
    w = meters_per_min_limit * math.sin(math.pi / 8.0)
    
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, 1., 0., 0., 0], -h)]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, -1., 0., 0., 0], -h)]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, 0., 1., 0., 0], -h)]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, 0., -1., 0., 0], -h)]

    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, 1., 1., 0., 0], -(w + h))]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, -1., 1., 0., 0], -(w + h))]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, -1., -1., 0., 0], -(w + h))]
    #ha.new_transition(p3, error).condition_list = [LinearConstraint([0, 0, 1., -1., 0., 0], -(w + h))]

    guard_1 = LinearConstraint([-1.0, 0.0, 0.0, 0.0, 0.0], 100)
    trans = ha.new_transition(p2, p3)
    trans.condition_list = [inv2_1, inv2_2, inv2_3, inv2_4, inv2_5, inv2_6, inv2_7, inv2_8]

    guard_2 = LinearConstraint([0.0, 0.0, 0.0, 0.0, -1.0, 0], -abortmintime)
    trans = ha.new_transition(p2, passive)
    trans.condition_list = [guard_2]

    guard_3 = LinearConstraint([0.0, 0.0, 0.0, 0.0, -1.0, 0], -abortmintime)
    trans = ha.new_transition(p3, passive)
    trans.condition_list = [guard_3]

    return ha

def define_init_states(ha):
    '''returns a list of (mode, HyperRectangle)'''
    # Variable ordering: [x, y, vx, vy, t]

    rv = []

    r = HyperRectangle([(-925.0, -875.0), (-425.0, -375.0), \
                        (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (1.0, 1.0)])
    rv.append((ha.modes['P2'], r))

    return rv

def define_settings(stepsize):
    'get the hylaa settings object'

    plot_settings = PlotSettings()
    #plot_settings.plot_mode = PlotSettings.PLOT_IMAGE if plot else PlotSettings.PLOT_NONE
    plot_settings.plot_mode = PlotSettings.PLOT_NONE

    plot_settings.max_shown_polys = None
    plot_settings.xdim = 0
    plot_settings.ydim = 1
    plot_settings.label.title = "Spacecraft Rendezvous (Abort)"
    plot_settings.filename = "rendezvous4_abort.png"
    plot_settings.plot_size = (10, 10)
    plot_settings.label.big(size=32)
    plot_settings.label.x_label = '$x$'
    plot_settings.label.y_label = '$y$'

    plot_settings.label.axes_limits = [-50, 50, -50, 50]
    #plot_settings.label.axes_limits = [-200, 200, -100, 100]

    y = 57.735
    lines1 = [(-100, y), (-100, -y), (0, 0), (-100, y)]
    lines2 = [(-5, -5), (5, -5), (5, 5), (-5, 5), (-5, -5)] 
    plot_settings.extra_lines = [lines1, lines2]

    #plot_settings.label.axes_limits = (-100, 100, -50, 0)

    settings = HylaaSettings(step=stepsize, max_time=250.0, plot_settings=plot_settings)

    settings.aggregation = False
    settings.deaggregation = False

    settings.skip_step_times = True
    settings.print_output = False

    settings.simulation.threads = 1

    return settings

def hylaa_v1_single(r, stepsize):
    'Runs hylaa with the given settings, returning the seconds elapsed'

    settings = define_settings(stepsize)

    ha = define_ha(r)
    init = define_init_states(ha)

    print "Running hylaa with r={} stepsize={}...".format(r, stepsize),

    engine = HylaaEngine(ha, settings)
    engine.run(init)

    rv = engine.result.time
    Timers.reset()

    print rv

    return rv

def hylaa_v1(r, stepsize=1.0):
    '''
    run hylaa with the given settings. will run multiple times if time is small
    '''

    print "small was 0.5"
    small = 0.5

    rv = hylaa_v1_single(r, stepsize=stepsize)

    if rv < small:
        # run four more times and return average of middle three runs
        times = [rv]

        for _ in range(4):
            rv = hylaa_v1_single(r, stepsize=stepsize)
            times.append(rv)

        times.sort()

        rv = sum(times[1:4]) / 3.0

        print "Runtimes: {}, returning {}".format(times, rv)

    return rv

def measure_hylaa_v1(timeout=60):
    'run hylaa measurements, saving to ../plot/data_hylaa_unagg.dat'

    print("Running hylaa_unagg measurements")
    start = time.time()
    step = 5
    params = [1, 0.5, 0.25]

    with open('../plot/data_hylaa_unagg.dat', 'w') as f:
        for stepsize in params:

            f.write('" step={}"\n'.format(stepsize))

            for rad in range(0, 75, step):
                sec = hylaa_v1(rad, stepsize=stepsize)

                f.write("{} {}\n".format(rad, sec))

                if sec > timeout: # don't need binary search here because we have a valid measurement
                    break

            f.write("\n\n")
            print("") # newline

    print "Finished hylaa_unagg in {} minutes".format((time.time() - start) / 60.0)

if __name__ == '__main__':
    measure_hylaa_v1(timeout=60)
