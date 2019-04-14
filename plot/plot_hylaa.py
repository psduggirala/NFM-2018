'''
Spacecraft Rendezvous System described in
"Verifying safety of an autonomous spacecraft rendezvous mission" by Nicole Chan and Sayan Mitra

This model was used in ARCHCOMP-18, specifically the variant with the abort maneuver.

Stanley Bak, July 2018
'''

import math

from matplotlib import collections

from hylaa.hybrid_automaton import HybridAutomaton
from hylaa.settings import HylaaSettings, PlotSettings, LabelSettings
from hylaa.core import Core
from hylaa.stateset import StateSet
from hylaa import lputil, aggstrat
from hylaa.aggstrat import Aggregated

def make_automaton(abortmin, abortmax):
    'make the hybrid automaton'

    safe = True
    ha = HybridAutomaton('Spacecraft Rendezvous with Abort')
    rad = 5.0

    passive_min_time = abortmin
    passive_max_time = abortmax
    
    ############## Modes ##############
    p2 = ha.new_mode('Far')
    a_mat = [\
        [0.0, 0.0, 1.0, 0.0, 0.0, 0], \
        [0.0, 0.0, 0.0, 1.0, 0.0, 0], \
        [-0.057599765881773, 0.000200959896519766, -2.89995083970656, 0.00877200894463775, 0.0, 0], \
        [-0.000174031357370456, -0.0665123984901026, -0.00875351105536225, -2.90300269286856, 0.0, 0.0], \
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], \
        [0, 0, 0, 0, 0, 0]]
    inv_mat = [[0.0, 0.0, 0.0, 0.0, 1.0, 0], [1.0, 0.0, 0.0, 0.0, 0.0, 0]]
    inv_rhs = [125.0, -100.0]
    p2.set_dynamics(a_mat)
    p2.set_invariant(inv_mat, inv_rhs)


    p3 = ha.new_mode('Approaching')
    a_mat = [\
        [0.0, 0.0, 1.0, 0.0, 0.0, 0], \
        [0.0, 0.0, 0.0, 1.0, 0.0, 0], \
        [-0.575999943070835, 0.000262486079431672, -19.2299795908647, 0.00876275931760007, 0.0, 0], \
        [-0.000262486080737868, -0.575999940191886, -0.00876276068239993, -19.2299765959399, 0.0, 0], \
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.],\
        [0, 0, 0, 0, 0, 0]]
    inv_mat = [\
        [-1.0, 0., 0., 0., 0., 0], \
        [1.0, 0., 0., 0., 0., 0], \
        [0., -1.0, 0., 0., 0., 0], \
        [0., 1.0, 0., 0., 0., 0], \
        [-1.0, -1.0, 0.0, 0.0, 0.0, 0], \
        [1.0, 1.0, 0.0, 0.0, 0.0, 0], \
        [-1.0, 1.0, 0., 0., 0., 0], \
        [1.0, -1.0, 0., 0., 0., 0], \
        [0., 0., 0., 0., 1., 0]]
    inv_rhs = [100, 100, 100, 100, 141.1, 141.1, 141.1, 141.1, passive_max_time]
    p3.set_dynamics(a_mat)
    p3.set_invariant(inv_mat, inv_rhs)


    passive = ha.new_mode('Abort')
    a_mat = [\
         [0, 0, 1, 0, 0, 0], \
         [0, 0, 0, 1, 0, 0], \
         [0.0000575894721132000, 0, 0, 0.00876276, 0, 0], \
         [0, 0, -0.00876276, 0, 0, 0], \
         [0, 0, 0, 0, 0, 1.], \
         [0, 0, 0, 0, 0, 0]]
    passive.set_dynamics(a_mat)
    
    error = ha.new_mode('Error')

    ############## Normal Transitions ##############
    t1 = ha.new_transition(p2, p3)
    guard_mat = [\
        [-1.0, 0., 0., 0., 0., 0], \
        [1.0, 0., 0., 0., 0., 0], \
        [0., -1.0, 0., 0., 0., 0], \
        [0., 1.0, 0., 0., 0., 0], \
        [-1.0, -1.0, 0.0, 0.0, 0.0, 0], \
        [1.0, 1.0, 0.0, 0.0, 0.0, 0], \
        [-1.0, 1.0, 0., 0., 0., 0], \
        [1.0, -1.0, 0., 0., 0., 0]]
    guard_rhs = [100, 100, 100, 100, 141.1, 141.1, 141.1, 141.1]
    t1.set_guard(guard_mat, guard_rhs)

    ha.new_transition(p2, passive).set_guard([[0.0, 0.0, 0.0, 0.0, -1.0, 0]], [-passive_min_time])

    ha.new_transition(p3, passive).set_guard([[0.0, 0.0, 0.0, 0.0, -1.0, 0]], [-passive_min_time])

    ############## Error Transitions ##############
    # In the aborting mode, the vehicle must avoid the target, which is modeled as a box B with
    # 0.2m edge length and the center placed as the origin

    # unsafe rad: 1.0
    #rad = 1.0
    
    t = ha.new_transition(passive, error)
    guard_mat = [ \
        [1, 0, 0., 0., 0., 0], \
        [-1, 0, 0., 0., 0., 0], \
        [0, 1., 0., 0., 0., 0], \
        [0, -1., 0., 0., 0., 0]]
    guard_rhs = [rad, rad, rad, rad]
    t.set_guard(guard_mat, guard_rhs)

    #In the rendezvous attempt the spacecraft must remain within the lineof-sight
    #cone L = {[x, y]^T | (x >= -100m) AND (y >= x*tan(30)) AND (-y >= x*tan(30))}
    #ha.new_transition(p3, error).set_guard([[1, 0, 0., 0., 0., 0]], [-100])
    #ha.new_transition(p3, error).set_guard([[-0.57735, 1, 0, 0., 0., 0]], [0])
    #ha.new_transition(p3, error).set_guard([[-0.57735, -1, 0., 0., 0., 0]], [0])

    # sqrt(vx^2 + vy^2) should stay below 0.055 m/SECOND (time in model is in MINUTES)
    # to make the model unsafe, try changing this to 0.05
    meters_per_sec_limit = 0.055 if safe else 0.05
    meters_per_min_limit = meters_per_sec_limit * 60
    h = meters_per_min_limit * math.cos(math.pi / 8.0)
    w = meters_per_min_limit * math.sin(math.pi / 8.0)
    
    #ha.new_transition(p3, error).set_guard([[0, 0, 1., 0., 0., 0]], [-h])
    #ha.new_transition(p3, error).set_guard([[0, 0, -1., 0., 0., 0]], [-h])
    #ha.new_transition(p3, error).set_guard([[0, 0, 0., 1., 0., 0]], [-h])
    #ha.new_transition(p3, error).set_guard([[0, 0, 0., -1., 0., 0]], [-h])

    #ha.new_transition(p3, error).set_guard([[0, 0, 1., 1., 0., 0]], [-(w + h)])
    #ha.new_transition(p3, error).set_guard([[0, 0, -1., 1., 0., 0]], [-(w + h)])
    #ha.new_transition(p3, error).set_guard([[0, 0, -1., -1., 0., 0]], [-(w + h)])
    #ha.new_transition(p3, error).set_guard([[0, 0, 1., -1., 0., 0]], [-(w + h)])
    
    return ha

def make_init(ha):
    'make the initial states'

    p2 = ha.modes['Far']
    box = [[-925.0, -875.0], [-425.0, -375.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0]]

    # move init y up 300
    #box[0][0] += 400
    #box[0][1] += 400
    
    init_lpi = lputil.from_box(box, p2)
    init_list = [StateSet(init_lpi, p2)]

    return init_list

def make_settings(stepsize, aggregation, axis_limits, filename):
    'make the reachability settings object'

    # see hylaa.settings for a list of reachability settings
    settings = HylaaSettings(stepsize, 250.0) # step: 0.1, bound: 200.0

    settings.plot.filename = filename

    #settings.approx_model = HylaaSettings.APPROX_LGG

    if aggregation == 'deagg':
        settings.stop_on_aggregated_error = False
        #settings.process_urgent_guards = True

        settings.aggstrat.deaggregate = True # use deaggregation
        settings.aggstrat.deagg_preference = Aggregated.DEAGG_LEAVES_FIRST
    elif aggregation == 'unagg':

        # dont use aggregation
        settings.aggstrat = aggstrat.Unaggregated()
    else:
        assert aggregation == 'agg', f"unknown aggregation value {aggregation}"

    settings.stdout = HylaaSettings.STDOUT_NONE

    #settings.plot.plot_mode = PlotSettings.PLOT_IMAGE
    #settings.plot.filename = "rendezvous_full_passivity.png"
    settings.plot.plot_size = (8, 4)

    settings.plot.video_pause_frames = 2
    settings.plot.plot_mode = PlotSettings.PLOT_IMAGE
    #settings.plot.filename = "rendezvous_full_passivity.mp4"

    # single plot
    settings.plot.xdim_dir = [0] * 1
    settings.plot.ydim_dir = [1] * 1
    settings.plot.grid = False
    settings.plot.label = []
    settings.plot.extra_collections = []

    ls = LabelSettings()
    settings.plot.label.append(ls)

    ls.big(size=24)
    ls.title = ''

    ls.x_label = '$x$'
    ls.y_label = '$y$'

    y = 57.735
    line = [(-100, y), (-100, -y), (0, 0), (-100, y)]
    c1 = collections.LineCollection([line], animated=True, colors=('gray'), linewidths=(1), linestyle='dashed')

    rad = 5.0
    line = [(-rad, -rad), (-rad, rad), (rad, rad), (rad, -rad), (-rad, -rad)]
    c2 = collections.LineCollection([line], animated=True, colors=('red'), linewidths=(2))

    settings.plot.extra_collections.append([c1, c2])

    settings.plot.label[0].axes_limits = axis_limits

    return settings

def plot_hylaa(abortmin, abortmax, stepsize=1.0, aggregation='deagg', axis_limits=None, filename='out.png'):
    'run hylaa a single time'

    ha = make_automaton(abortmin, abortmax)

    init_states = make_init(ha)

    settings = make_settings(stepsize, aggregation, axis_limits, filename)

    print(f"running {aggregation} with step {stepsize}")
    result = Core(ha, settings).run(init_states)

    safe = not result.has_concrete_error

    if aggregation == 'agg':
        safe = not result.has_aggregated_error
    
    secs = result.top_level_timer.total_secs

    print(f"finished in {secs} secs")

    return [safe, secs]

if __name__ == "__main__":
    # 5 secs
    #plot_hylaa(140, 140, stepsize=1.0, aggregation='unagg', axis_limits=[-20, 20, -15, 15], filename='step_a_1.0.png')

    # 15 secs
    #plot_hylaa(140, 140, stepsize=0.5, aggregation='unagg', axis_limits=[-20, 20, -15, 15], filename='step_b_0.5.png')

    # 50 secs
    #plot_hylaa(140, 140, stepsize=0.25, aggregation='unagg', axis_limits=[-20, 20, -15, 15], filename='step_c_0.25.png')

    # 50 secs
    plot_hylaa(140, 140, stepsize=0.25, aggregation='unagg', axis_limits=[-1000, 200, -500, 100], filename='step_d_0.25.png')

    # 5 secs
    #plot_hylaa(140, 140, stepsize=1.0, aggregation='deagg', axis_limits=[-20, 20, -15, 15], filename='deagg_1.0.png')
