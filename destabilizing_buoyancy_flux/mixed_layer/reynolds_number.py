import os, sys; sys.path.append(os.getenv('DEDALES', os.path.join("..", "..")))

import time, logging
import numpy as np

from mpi4py import MPI
from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

logging.basicConfig(datefmt="%I:%M:%S")
logger = logging.getLogger(__name__)

# Some convenient constants
second = 1.0
minute = 60*second
hour   = 60*minute
day    = 24*hour

# Physical constants
α  = 2.0e-4   # Thermal expansion coefficient [K⁻¹]
β  = 8.0e-4   # Saline expansion coefficient [no units]
g  = 9.81     # Graviational acceleration [m s⁻²]
ρ0 = 1028.1   # Reference density [kg m⁻³]
cP = 3993.0   # Specific heat of oceanic water [J kg⁻¹ K⁻¹]

def float_2_nice_str(a):
    before = "{:.0f}".format(a)
    after = "{:.15f}".format(a-np.floor(a))
    return "{}p{}".format(before, after[2:-1].rstrip('0'))

def identifier(model, closure=None): 
    if closure is None: closure_name = 'DNS'
    else:               closure_name = closure.__class__.__name__
    return "case{}_{:s}".format(model.case, closure_name)

# Parameters for investigating Reynolds number effects on convection
#
#
#
parameters = {
    '0': {'Lx': 50.0, 'Lz': 25.0, 'Fb': 1e-8, 'N0': 1/1000.0, 'ν': 1.05e-6, 'κ': 1.43e-7},  
    '1': {'Lx': 50.0, 'Lz': 25.0, 'Fb': 1e-8, 'N0': 1/1000.0, 'ν': 1.00e-2, 'κ': 1.00e-2},
    '2': {'Lx': 50.0, 'Lz': 25.0, 'Fb': 1e-8, 'N0': 1/1000.0, 'ν': 5.00e-3, 'κ': 5.00e-3},
    '3': {'Lx': 50.0, 'Lz': 25.0, 'Fb': 1e-8, 'N0': 1/1000.0, 'ν': 2.00e-3, 'κ': 2.00e-3},
    '4': {'Lx': 50.0, 'Lz': 25.0, 'Fb': 1e-8, 'N0': 1/1000.0, 'ν': 1.00e-3, 'κ': 1.00e-3},
}

if sys.argv[-1] == 'debug':
    debug = True
else:
    debug = False

try:
    case = sys.argv[2]
    params = parameters[case]
except:
    case = '1'
    params = parameters[case]
    logger.info("Case not found! Using case {}.".format(case))

try:
    closure_name = sys.argv[1]
    if closure_name is 'DNS':
        closure = None
    elif closure_name == 'ConstantSmagorinsky':
        closure = getattr(dedaLES, closure_name)()
        dt_safety = 0.5
    else:
        closure = getattr(dedaLES, closure_name)()
        dt_safety = 1.0
except:
    closure_name = 'DNS'
    closure = None
    logger.info("Closure not found! Using {}.".format(closure_name))
    dt_safety = 0.5

# Main parameters
nx = ny = 128 # x,y resolution 
nz = 64 
Lx = Ly = params['Lx'] # x,y extent [m]
Lz = params['Lz']

surface_buoyancy_flux = params['Fb']    # Buoyancy flux into ocean [m² s⁻³]
initial_N = params['N0']                # Initial buoyancy frequency [s⁻¹]
ν = params['ν']
κ = params['κ']
initial_h = 10.0                        # Initial mixed layer depth [m]

# Physical parameters
turb_vel_scale          = (Lz*surface_buoyancy_flux)**(1/3)         # Domain turbulent velocity scale [m s⁻¹]
noise_amplitude         = 0.1*turb_vel_scale                        # Noise amplitude [m s⁻¹]

surface_bz              = -surface_buoyancy_flux/κ                  # [s⁻²]
initial_dt              = 1e-4 / np.sqrt(-surface_bz)

surface_heating         = -surface_buoyancy_flux*ρ0*cP/(α*g)        # [W m⁻²]
kolmogorov_length_scale = (ν**3/surface_buoyancy_flux)**(1/4)       # Kolmogorov length scale
initial_N2              = initial_N**2                              # Initial buoyancy gradient [s⁻²]
erosion_time_scale      = initial_N2*Lz**2/surface_buoyancy_flux    # Time-scale for stratification erosion

# Numerical parameters
dt_cadence       = 10
stats_cadence    = 100
averages_cadence = 10
analysis_cadence = 100
run_time         = day
max_writes       = 100

if debug:
    nx = ny = nz = 8
    dt_cadence = np.inf
    #initial_dt = 1e-16
    run_time = 10*initial_dt
    stats_cadence = analysis_cadence = averages_cadence = 1

# Construct model
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, closure=closure, case=int(case),
                                      surface_buoyancy_flux=surface_buoyancy_flux, surface_bz=surface_bz, initial_N2=initial_N2)

Δx = Lx/nx
Δy = Ly/ny
Δz_min, Δz_max = dedaLES.grid_stats(model, 2)
Δmin = min(Δx, Δy, Δz_min)

logger.info("""\n
    *** Convection into a linearly stratified fluid ***

                       Simulation info
                       ---------------

              surface heating : {:.2e} W m⁻²
        surface buoyancy flux : {:.2e} m² s⁻³
                          1/N : {:.2e} s
                   initial dt : {:.2e} s
                     run time : {:.2e} s

                           Lx : {:.1f} m
                           Ly : {:.1f} m
                           Lz : {:.1f} m

                           nx : {:d}
                           ny : {:d}
                           nz : {:d}

            turbulent w-scale : {:.2e} m s⁻¹
           erosion time-scale : {:.2e} s
             Kolmogorov scale : {:.2e} m
                    x-spacing : {:.2e} m 
                    y-spacing : {:.2e} m
           z-spacing min, max : {:.2e} m, {:.2e} m
                      closure : {}

    """.format(surface_heating, surface_buoyancy_flux, 1/initial_N, initial_dt, run_time, 
                Lx, Ly, Lz, nx, ny, nz,
                turb_vel_scale, erosion_time_scale, kolmogorov_length_scale, Δx, Δy, Δz_min, Δz_max,
                closure_name)
)

# Boundary conditions
model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="surface_bz") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="initial_N2")
model.build_solver(timestepper='RK222')

# Initial condition
noise = noise_amplitude * dedaLES.random_noise(model.domain) * model.z * (Lz - model.z) / Lz**2
initial_b = initial_N2*model.z
model.set_fields(
    u = noise,
    v = noise,
    b = initial_b + np.sqrt(initial_N2) * noise
)

model.stop_at(sim_time=run_time)

dt_gizmo = dedaLES.TimeStepGizmo(model.solver, initial_dt=initial_dt, cadence=dt_cadence, max_change=1.5, safety=dt_safety)
dt_gizmo.add_velocities(('u', 'v', 'w'))
dt_gizmo.add_diffusivity('ν_sgs')
dt_gizmo.add_diffusivity('κb_sgs')

stats = flow_tools.GlobalFlowProperty(model.solver, cadence=stats_cadence)

stats.add_property("w*b", name="wb")
stats.add_property("ε + ε_sgs", name="epsilon")
stats.add_property("ν_sgs", name="eddy viscosity")
stats.add_property("κb_sgs", name="eddy diffusivity")
stats.add_property("χ + χ_sgs", name="chi")
stats.add_property("w*w", name="w_sq")
stats.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')

analysis = model.solver.evaluator.add_file_handler("canonical_analysis_{}".format(identifier(model, closure=closure)), 
                                                   iter=analysis_cadence, max_writes=max_writes)
analysis.add_system(model.solver.state, layout='g')
analysis.add_task("interp(b, y=0)", scales=1, name='b midplane')
analysis.add_task("interp(u, y=0)", scales=1, name='u midplane')
analysis.add_task("interp(v, y=0)", scales=1, name='v midplane')
analysis.add_task("interp(w, y=0)", scales=1, name='w midplane')

averages = model.solver.evaluator.add_file_handler("canonical_averages_{}".format(identifier(model, closure=closure)), 
                                                   iter=averages_cadence, max_writes=max_writes)
averages.add_task("integ(integ(u, 'x'), 'y')", scales=1, name='avg u')
averages.add_task("integ(integ(v, 'x'), 'y')", scales=1, name='avg v')
averages.add_task("integ(integ(w, 'x'), 'y')", scales=1, name='avg w')
averages.add_task("integ(integ(b, 'x'), 'y')", scales=1, name='avg b')

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    log_time = time.time()  

    while model.solver.ok:
        dt = dt_gizmo.compute_dt()
        model.solver.step(dt)
        if model.time_to_log(stats_cadence): 
            compute_time = time.time() - log_time
            log_time = time.time()

            dt_sgs = Δmin**2 / stats.max("eddy viscosity")

            logger.info("""i: {:d}, t: {:.3f} hr ({:.1f} %), twall: {:.1f} s, dt: {:.2e} s, max Re {:.0f}, max sqrt(w^2): {:.2e}, 
    max ε: {:.2e}, <ε>: {:.2e}, <wb>: {:.2e}, <χ>: {:.2e}, max ν_sgs: {:.2e}, max κ_sgs: {:.2e}, dt_sgs: {:.2e}""".format( 
                model.solver.iteration, model.solver.sim_time/hour, 0.01*model.solver.sim_time/run_time, compute_time, dt, stats.max("Re"),
                np.sqrt(stats.max("w_sq")), stats.max("epsilon"), stats.volume_average("epsilon"),
                stats.volume_average("wb"), stats.volume_average("chi"), stats.max("eddy viscosity"), stats.max("eddy diffusivity"), dt_sgs
            ))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise

finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %model.solver.iteration)
    logger.info('Sim end time: %f' %model.solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * model.domain.dist.comm_cart.size))
