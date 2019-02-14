import os, sys; sys.path.append(os.getenv('DEDALES', os.path.join("..", "..")))

import time, logging
import numpy as np

from mpi4py import MPI
from numpy import pi
from dedalus.extras import flow_tools

import dedaLES

logging.basicConfig(datefmt="%I:%M:%S")
logger = logging.getLogger(__name__)

# Setup
debug = False

if len(sys.argv) is 1:
    closure = None
    closure_name = 'DNS'
else:
    closure_name = sys.argv[1]
    try:
        if closure_name is 'DNS':
            closure = None
        else:
            closure = getattr(dedaLES, sys.argv[1])()
    except:
        logger.info("Closure '{}' not found! Running in debug mode.".format(closure_name))
        closure = None
        debug = True

# Some convenient constants
second = 1.0
minute = 60*second
hour   = 60*minute
day    = 24*hour

# Physical constants
κ  = 1.43e-7  # Thermal diffusivity [m² s⁻¹]
ν  = 1.05e-6  # Viscosity [m² s⁻¹]
α  = 2.0e-4   # Thermal expansion coefficient [K⁻¹]
β  = 8.0e-4   # Thermal expansion coefficient [K⁻¹]
g  = 9.81     # Graviational acceleration [m s⁻²]
ρ0 = 1028.1   # Reference density [kg m⁻³]
cP = 3993.0   # Specific heat of oceanic water [J kg⁻¹ K⁻¹]

def heat_flux(buoyancy_flux):
    return ρ0*cP/(α*g) * buoyancy_flux

def buoyancy_flux(heat_flux):
    return α*g/(ρ0*cP) * heat_flux

def float_2_nice_str(a):
    before = "{:.0f}".format(a)
    after = "{:.15f}".format(a-np.floor(a))
    return "{}p{}".format(before, after[2:-1].rstrip('0'))

def identifier(model, closure=None): 
    if closure is None: closure_name = 'DNS'
    else:               closure_name = closure.__class__.__name__
    return "nx{:d}_ny{:d}_nz{:d}_F{}_Ninv{:.0f}_{:s}".format(
            model.nx, model.ny, model.nz, float_2_nice_str(-model.surface_buoyancy_flux), 1/np.sqrt(initial_N2), closure_name)

# Main parameters
nx = ny = nz = 128   # x,y resolution 
Lx = Ly = Lz = 64.0  # x,y,z extent [m]

surface_buoyancy_flux = -1e-9  # Buoyancy flux into ocean [m² s⁻³]
initial_N = 1/500.0    # Initial buoyancy frequency [s⁻¹]
initial_h = 10.0       # Initial mixed layer depth [m]

# Physical parameters
turb_vel_scale          = (-Lz*surface_buoyancy_flux)**(1/3)        # Domain turbulent velocity scale [m s⁻¹]
noise_amplitude         = 0.01*turb_vel_scale                       # Noise amplitude [m s⁻¹]

surface_bz              = surface_buoyancy_flux/κ                   # [s⁻²]
initial_dt              = 1e-2 / np.sqrt(-surface_bz)

surface_heat_flux       = surface_buoyancy_flux*ρ0*cP/(α*g)         # [W m⁻²]
kolmogorov_length_scale = (-ν**3/surface_buoyancy_flux)**(1/4)      # Kolmogorov length scale
initial_N2              = initial_N**2                              # Initial buoyancy gradient [s⁻²]
erosion_time_scale      = -initial_N2*Lz**2/surface_buoyancy_flux   # Time-scale for stratification erosion

# Numerical parameters
CFL_cadence      = 10
stats_cadence    = 100
averages_cadence = 10
analysis_cadence = 100
run_time         = 2*hour
max_writes       = 1000

if debug:
    nx = ny = nz = 8
    CFL_cadence = np.inf
    initial_dt = 1e-16
    run_time = 10*initial_dt
    stats_cadence = analysis_cadence = averages_cadence = 1

# Construct model
model = dedaLES.BoussinesqChannelFlow(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz, ν=ν, κ=κ, closure=closure,
                                      surface_buoyancy_flux=surface_buoyancy_flux, surface_bz=surface_bz, initial_N2=initial_N2)

Δx = Lx/nx
Δy = Ly/ny
Δz_min, Δz_max = dedaLES.grid_stats(model, 2)

logger.info("""\n
    *** Convection into a linearly stratified fluid ***

                       Simulation info
                       ---------------

            surface heat_flux : {:.2e} W m⁻²
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

    """.format(surface_heat_flux, surface_buoyancy_flux, 1/initial_N, initial_dt, run_time, 
               Lx, Ly, Lz, nx, ny, nz,
               turb_vel_scale, erosion_time_scale, kolmogorov_length_scale, Δx, Δy, Δz_min, Δz_max)
)

# Boundary conditions
model.set_bc("no penetration", "top", "bottom")
model.set_bc("free slip", "top", "bottom")
model.set_tracer_gradient_bc("b", "top", gradient="surface_bz") #*tanh(t/tb)")
model.set_tracer_gradient_bc("b", "bottom", gradient="initial_N2")
model.build_solver(timestepper='SBDF3')

# Initial condition
noise = noise_amplitude * dedaLES.random_noise(model.domain) * model.z * (Lz - model.z) / Lz**2
initial_b = initial_N2*model.z
model.set_fields(
    u = noise,
    v = noise,
    b = initial_b + np.sqrt(initial_N2) * noise
    #w = noise,
)

model.stop_at(sim_time=run_time)

CFL = flow_tools.CFL(model.solver, initial_dt=initial_dt, cadence=CFL_cadence, max_change=1.5, safety=0.5)
CFL.add_velocities(('u', 'v', 'w'))

stats = flow_tools.GlobalFlowProperty(model.solver, cadence=stats_cadence)

stats.add_property("w*b", name="wb")
stats.add_property("ε + ε_sgs", name="epsilon")
stats.add_property("χ + χ_sgs", name="chi")
stats.add_property("w*w", name="w_sq")
stats.add_property("sqrt(u*u + v*v + w*w) / ν", name='Re')

analysis = model.solver.evaluator.add_file_handler("mixed_layer_analysis_{}".format(identifier(model, closure=closure)), 
                                                   iter=analysis_cadence, max_writes=max_writes)
analysis.add_system(model.solver.state, layout='g')
analysis.add_task("interp(b, y=0)", scales=1, name='b midplane')
analysis.add_task("interp(u, y=0)", scales=1, name='u midplane')
analysis.add_task("interp(v, y=0)", scales=1, name='v midplane')
analysis.add_task("interp(w, y=0)", scales=1, name='w midplane')

averages = model.solver.evaluator.add_file_handler("mixed_layer_averages_{}".format(identifier(model, closure=closure)), 
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
        dt = CFL.compute_dt()
        model.solver.step(dt)
        if model.time_to_log(stats_cadence): 
            compute_time = time.time() - log_time
            log_time = time.time()

            logger.info("""i: {:d}, t: {:.3f} hr, twall: {:.1f} s, dt: {:.2f} s, max Re {:.0f} 
    max sqrt(w^2): {:.2e}, max ε: {:.2e}, <ε>: {:.2e}, <wb>: {:.2e}, <χ>: {:.2e}""".format( 
        model.solver.iteration, model.solver.sim_time/hour, compute_time, dt, stats.max("Re"),
        np.sqrt(stats.max("w_sq")), stats.max("epsilon"), stats.volume_average("epsilon"),
        stats.volume_average("wb"), stats.volume_average("chi")
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
