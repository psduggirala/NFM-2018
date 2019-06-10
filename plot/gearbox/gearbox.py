'''
Gearbox benchmark from:
"Motor-Transmission Drive System: a Benchmark Example for Safety Verification"
by Hongxu Chen, Sayan Mitra and Guangyu Tian in ARCH 2015

This model was also used in ARCHCOMP-17 and 18.

This model shows how to use hylaa.symbolic to construct the dynamics and reset / guard conditions from 
string expressions.

Unique features:
- Performs both reachability analysis and simulation on the model
- Uses deaggregation to find worst-case meshing time and cumulative impulse force
- Concurrently creates two plots at once (x/y and impuse force/time)

Stanley Bak, Nov 2018
'''

import math

import numpy as np

from matplotlib import collections

from hylaa.hybrid_automaton import HybridAutomaton
from hylaa.settings import HylaaSettings, PlotSettings, LabelSettings
from hylaa.core import Core
from hylaa.stateset import StateSet
from hylaa import lputil, aggstrat, symbolic
from hylaa.aggstrat import Aggregated

from gearbox_sim import get_worst_case_sim

def make_automaton(theta_deg, maxi=20):
    'make the hybrid automaton'

    ha = HybridAutomaton('Gearbox')

    # variables
    variables = ["px", "py", "vx", "vy", "I"]
    derivatives = ["vx", "vy", "Fs/ms", "-Rs*Tf/Jg2", "0"]

    constant_dict = {
        "zeta": 0.9,
        "ms": 3.4,
        "mg2": 18.1,
        "Jg2": 0.7,
        "Rs": 0.08,
        "theta": theta_deg * math.pi / 180, #0.628318530717959,
        "deltap": -0.003,
        "Fs": 70,
        "Tf": 1
        }

    ############## Modes ##############
    move_free = ha.new_mode('move_free')
    #
    a_mat = symbolic.make_dynamics_mat(variables, derivatives, constant_dict, has_affine_variable=True)
    move_free.set_dynamics(a_mat)

    invariant = f"px<=deltap & py<=-px*0.726542528005361 & py>=px*0.726542528005361 & I <= {maxi}"
    mat, rhs = symbolic.make_condition(variables, invariant.split('&'), constant_dict, has_affine_variable=True)
    move_free.set_invariant(mat, rhs)

    meshed = ha.new_mode('meshed')
    a_mat = np.zeros((6, 6))
    meshed.set_dynamics(a_mat)

    # error mode
    error = ha.new_mode('error')

    ############## Cyclic Transitions ##############
    # t1
    t1 = ha.new_transition(move_free, move_free, "t1")
    t1_guard = "py>=-px*0.726542528005361 & vx*0.587785252292473+vy*0.809016994374947>=0 & px <= deltap"
    mat, rhs = symbolic.make_condition(variables, t1_guard.split('&'), constant_dict, has_affine_variable=True)
    t1.set_guard(mat, rhs)

    # projection of a point onto a plane
    # if <x, y> is the point and <nx, ny> is the normal to the plane with mag 1
    # proj = <x, y> - (<x, y> dot <nx, ny>) * <nx, ny>
    nx = math.cos(math.pi / 2 - constant_dict['theta'])
    ny = math.sin(math.pi / 2 - constant_dict['theta'])

    i_reset = "I+(vx*0.587785252292473+vy*0.809016994374947)*(zeta+1)*ms*mg2/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))" 
    resets = [f"px - (px*{nx} + py*{ny}) * {nx}",
              f"py - (px*{nx} + py*{ny}) * {ny}",
              "(vx*(ms*0.809016994374947*0.809016994374947-mg2*zeta*0.587785252292473*0.587785252292473)+vy*(-(zeta+1)*mg2*0.587785252292473*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))", 
              "(vx*(-(zeta+1)*ms*0.587785252292473*0.809016994374947)+vy*(mg2*0.587785252292473*0.587785252292473-ms*zeta*0.809016994374947*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))", 
              i_reset]
    reset_mat = symbolic.make_reset_mat(variables, resets, constant_dict, has_affine_variable=True)
    t1.set_reset(reset_mat)

    # T1 to error mode
    t1_to_error = ha.new_transition(move_free, error, "t1_to_error")
    guard = f"{t1_guard} & {i_reset}>={maxi}"
    mat, rhs = symbolic.make_condition(variables, guard.split('&'), constant_dict, has_affine_variable=True)
    t1_to_error.set_guard(mat, rhs)
    t1_to_error.set_reset(reset_mat)

    # t2
    t2 = ha.new_transition(move_free, move_free, "t2")
    t2_guard = "py<=px*0.726542528005361 & vx*0.587785252292473-vy*0.809016994374947>=0 & px <= deltap"
    mat, rhs = symbolic.make_condition(variables, t2_guard.split('&'), constant_dict, has_affine_variable=True)
    t2.set_guard(mat, rhs)

    ny = -ny

    i_reset = "I+(vx*0.587785252292473-vy*0.809016994374947)*(zeta+1)*ms*mg2/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))"
    resets = [f"px - (px*{nx} + py*{ny}) * {nx}",
              f"py - (px*{nx} + py*{ny}) * {ny}",
              "(vx*(ms*0.809016994374947*0.809016994374947-mg2*zeta*0.587785252292473*0.587785252292473)+vy*((zeta+1)*mg2*0.587785252292473*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))",
              "(vx*((zeta+1)*ms*0.587785252292473*0.809016994374947)+vy*(mg2*0.587785252292473*0.587785252292473-ms*zeta*0.809016994374947*0.809016994374947))/(ms*(0.809016994374947*0.809016994374947)+mg2*(0.587785252292473*0.587785252292473))",
              i_reset]
    reset_mat = symbolic.make_reset_mat(variables, resets, constant_dict, has_affine_variable=True)
    t2.set_reset(reset_mat)

    # T2 to error mode
    t2_to_error = ha.new_transition(move_free, error, "t2_to_error")
    guard = f"{t2_guard} & {i_reset}>={maxi}"
    mat, rhs = symbolic.make_condition(variables, guard.split('&'), constant_dict, has_affine_variable=True)
    t2_to_error.set_guard(mat, rhs)
    t2_to_error.set_reset(reset_mat)

    #### transitions to meshed

    guards = ["px >= deltap & vx >= 0 & vy >= 0",
              "px >= deltap & vx >= 0 & vy <= 0",
              "px >= deltap & vx <= 0 & vy >= 0",
              "px >= deltap & vx <= 0 & vy <= 0"]

    i_resets = ["I+ms*vx+ms*vy", "I+ms*vx-ms*vy", "I-ms*vx+ms*vy", "I-ms*vx-ms*vy"]

    for i, (guard, i_reset) in enumerate(zip(guards, i_resets)):
        t3 = ha.new_transition(move_free, meshed, f"meshed_{i}")

        t3_guard = guard + f" & I <= {maxi}"
        mat, rhs = symbolic.make_condition(variables, t3_guard.split('&'), constant_dict, has_affine_variable=True)
        t3.set_guard(mat, rhs)

        resets = [f"px",
                  f"py",
                  "0",
                  "0",
                  i_reset]
        reset_mat = symbolic.make_reset_mat(variables, resets, constant_dict, has_affine_variable=True)
        t3.set_reset(reset_mat)

        # T3 to error mode
        t3_to_error = ha.new_transition(move_free, error, f"t3_to_error_{i}")
        t3_err_guard = guard + f" & I >= {maxi}"
        mat, rhs = symbolic.make_condition(variables, t3_err_guard.split('&'), constant_dict, has_affine_variable=True)
        t3_to_error.set_guard(mat, rhs)
            
    return ha

def make_init(ha, box):
    'make the initial states'

    mode = ha.modes['move_free']
    # px==-0.0165 & py==0.003 & vx==0 & vy==0 & I==0 & affine==1.0
    init_lpi = lputil.from_box(box, mode)

    
    #init_lpi = lputil.from_box([(-0.02, -0.02), (-0.005, -0.003), (0, 0), (0, 0), (0, 0), (1.0, 1.0)], mode)
    #start = [-0.02, -0.004213714568273684, 0.0, 0.0, 0.0, 1.0]
    #tol = 1e-7
    #init_lpi = lputil.from_box([(x - tol, x + tol) if i < 2 else (x, x) for i, x in enumerate(start)], mode)

    init_list = [StateSet(init_lpi, mode)]

    # why does 0.003-0.005 reach an error with i=30 for roots first but not leaves first?
 
    return init_list

def make_settings(theta_deg, box, tmax):
    'make the reachability settings object'

    # see hylaa.settings for a list of reachability settings
    settings = HylaaSettings(0.001, tmax) # step: 0.001, bound: 0.35

    #settings.stop_on_aggregated_error = False

    settings.skip_zero_dynamics_modes = False
    
    #settings.aggstrat = MyAggergated()
    settings.aggstrat.deaggregate = True # use deaggregation
    settings.aggstrat.deagg_preference = Aggregated.DEAGG_LEAVES_FIRST

    settings.stdout = HylaaSettings.STDOUT_VERBOSE

    #settings.plot.plot_mode = PlotSettings.PLOT_IMAGE #INTERACTIVE
    #settings.plot.filename = "gearbox.png"
    settings.plot.plot_mode = PlotSettings.PLOT_NONE
    settings.plot.plot_size = (8, 9)

    settings.plot.xdim_dir = 0
    settings.plot.ydim_dir = 1

    settings.plot.xdim_dir = [0, None]
    settings.plot.ydim_dir = [1, 4]
    settings.plot.label = [LabelSettings(), LabelSettings()]

    settings.plot.label[0].title = "Gearbox"
    settings.plot.label[0].title_size = 22
    
    #settings.plot.label.axes_limits = [-0.02, 0, -0.008, 0.004]
    settings.plot.label[0].axes_limits = [-0.02, 0, -0.01, 0.01]
    settings.plot.label[0].x_label = "X Position [m]"
    settings.plot.label[0].y_label = "Y Position [m]"
    settings.plot.label[0].label_size = 20

    settings.plot.label[1].axes_limits = [0.0, 0.35, 0.0, 40.0]
    settings.plot.label[1].x_label = "Time [s]"
    settings.plot.label[1].y_label = "Impact Impulse [N m]"
    settings.plot.label[1].label_size = 20

    lines = []

    # add coords
    theta = theta_deg * math.pi / 180
    deltap = 0.003
    b = 0.01
    
    y = deltap * math.tan(theta)
    x = b / (2 * math.tan(theta))

    p1 = (-x - deltap, y + b / 2)
    p2 = (-deltap, y)
    p3 = (-deltap, -y)
    p4 = (-x - deltap, -y - b / 2)
    
    lines.append([p1, p2, p3, p4])

    line_col = collections.LineCollection(lines, animated=True, colors=('gray'), linewidths=(2), linestyle='dashed')

    verts = [[box[0][0], box[1][0]], [box[0][1], box[1][0]], [box[0][1], box[1][1]], [box[0][0], box[1][1]]]

    poly_col = collections.PolyCollection([verts], animated=True, edgecolor=None, facecolor=(0., 1.0, 0., 0.5))

    settings.plot.extra_collections = [[line_col, poly_col], []]

    return settings

def run_hylaa(is_sim=False, tmax=0.1):
    'main entry point'

    theta_deg = 36
    maxi = 60

    ha = make_automaton(theta_deg, maxi)

    box = [(-0.017, -0.016), (-0.005, 0.005), (0, 0), (0, 0), (0, 0), (1.0, 1.0)]

    settings = make_settings(theta_deg, box, tmax)

    if not is_sim:
        # add the worst-case sim on top
        plot1_sim, plot2_sim = get_worst_case_sim()

        plot1_lc = collections.LineCollection([plot1_sim], animated=True, colors=('red'), linewidths=(1), zorder=6)
        settings.plot.extra_collections[0].append(plot1_lc)

        box_pts = [(box[0][0], box[1][0]), (box[0][1], box[1][0]), (box[0][1], box[1][1]), (box[0][0], box[1][1]),
                   (box[0][0], box[1][0])]
        plot1_box = collections.LineCollection([box_pts], animated=True, colors=('black'), linewidths=2, zorder=6)
        settings.plot.extra_collections[0].append(plot1_box)

        plot2_lc = collections.LineCollection([plot2_sim], animated=True, colors=('red'), linewidths=(1), zorder=6)
        settings.plot.extra_collections[1].append(plot2_lc)

        settings.plot.filename = "gearbox.png"
        settings.plot.plot_mode = PlotSettings.PLOT_IMAGE

        # settings for tiem measurements
        #settings.plot.plot_mode = PlotSettings.PLOT_NONE
        #settings.stdout = HylaaSettings.STDOUT_NORMAL
        #settings.skip_zero_dynamics_modes = True
        
        result = Core(ha, settings).run(make_init(ha, box))

        if result.counterexample:
            print(f"counterexample start: {result.counterexample[0].start}")
    else:
        init_mode = ha.modes['move_free']
        settings.aggstrat.sim_avoid_modes.append('meshed')
        settings.plot.sim_line_width = 1.0
        settings.plot.filename = "gearbox_sim.png"
        #settings.plot.filename = "gearbox_sim.mp4"

        settings.plot.plot_mode = PlotSettings.PLOT_IMAGE

        settings.random_seed = 0 # 1 for 100 sims no worst-case
        result = Core(ha, settings).simulate(init_mode, box, [84, 85]) # was 100
        #[84, 85] [18. 19]

        # plot -> segment -> sim_id -> list_of_points
        #print(f"sim_lines len: {len(result.sim_lines)}")
        #print(f"sim_lines[0] len: {len(result.sim_lines[0])}")
        #print(f"sim_lines[0][0] len: {len(result.sim_lines[0][0])}\n")
        #print(f"sim_lines[0][0][0] len: {len(result.sim_lines[0][0][0])}\n")

        # this code prints out one counterexample simulation to the terminal
        if False:
            sim_lines = result.sim_lines

            for plot in [0, 1]:
                for segment in sim_lines[plot]:
                    sim_id = 0

                    for line in segment[sim_id]:
                        print(f"pts.append({line})")

                print("\n")

        #print(f"sim_lines[0][0][0][0]: {result.sim_lines[0][0][0][0]}")
        #print(f"sim_lines[0][0][0][1]: {result.sim_lines[0][0][0][1]}")
        #print(f"sim_lines[1][0][0][0]: {result.sim_lines[1][0][0][0]}")

if __name__ == "__main__":
    tmax = 0.35 # 0.35
    
    run_hylaa(is_sim=False, tmax=tmax)
    #run_hylaa(is_sim=True, tmax=tmax)
