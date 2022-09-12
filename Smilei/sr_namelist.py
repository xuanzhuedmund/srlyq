from math import pi, sqrt

l0 = 2. * pi             # laser wavelength [in code units]
t0 = l0                  # optical cycle
Lsim = [60*l0, 35.*l0]  # length of the simulation
Tsim = 55.*t0            # duration of the simulation
cpx = 152
cpy = 152                 #proportion about nx,ny
lambdar = 0.8e-6
c = 299792458
wr = 2*math.pi*c/lambdar
Main(
    geometry="2Dcartesian",

    interpolation_order=2,

    cell_length=[l0 / cpx , l0 / cpy],
    grid_length=Lsim,

    number_of_patches=[32, 8],  #the number of patches in each direction
                               # Each integer must be a power of 2 and
                               # the total number of patches must be
                               # greater than the number of MPI processes

    timestep=0.8*t0/(sqrt(cpx**2+cpy**2)),
    #timestep_over_CFL=0.8, #CFL=t0/sqrt(cpx**2+cpy**2)
    simulation_time=Tsim,

    EM_boundary_conditions=[
        ['silver-muller'],
        ['silver-muller'],
    ],

    random_seed=smilei_mpi_rank,
    reference_angular_frequency_SI=wr,
)

omega = 1.
a0 = 118*math.sqrt(2) #a0 of CP
waist = 1.5*l0 #The waist value. Transverse coordinate at which the field is at 1/e of its maximum value.
Zr = omega * waist**2/2. #x_R = pi * w_0^2 / lambda0, confocal parameter f(Rayleigh length)
focus = [10.*l0, 17.5*l0] #The X and Y positions of the laser focus.
polarization_phi = 0
ellipticity = 1
w  = math.sqrt(1./(1.+(focus[0]/Zr)**2))
invWaist2 = (w/waist)**2
coeff = -omega * focus[0] * w**2 / (2.*Zr**2)

dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
amplitudeY *= a0 * omega
amplitudeZ *= a0 * omega
delay_phase = [dephasing, 0]
def By(y,t):
    B = amplitudeY * w * math.exp( -invWaist2*(y-focus[1])**2 )**1 \
        * math.sin(omega*t - coeff*(y-focus[1])**2 + delay_phase[0])
    if t < t0: return t/t0 * B
    elif t< 9*t0: return B
    elif t< 10*t0: return (10 - t/t0) * B
    else: return 0
def Bz(y,t):
    B = amplitudeZ * w * math.exp( -invWaist2*(y-focus[1])**2 )**1 \
        * math.sin(omega*t - coeff*(y-focus[1])**2 + delay_phase[1])
    if t < t0: return t/t0 * B
    elif t< 9*t0: return B
    elif t< 10*t0: return (10 - t/t0) * B
    else: return 0

Laser(
    box_side           = "xmin",
    space_time_profile = [By, Bz]
)
# LaserGaussian2D(
#     box_side         = "xmin",
#     a0               = 69.08*math.sqrt(2),
#     omega            = 1.,
#     focus            = [10.*l0, 17.*l0],
#     waist            = 2.0*l0,
#     incidence_angle  = 0.,
#     polarization_phi = 0.,
#     ellipticity      = 1.,
#     time_envelope    = ttrapezoidal(plateau=8*t0, slope1=t0, slope2=t0)
# )
n_C = 59.31
n_H = 19.77
target_min = 10*l0
target_max = 10.125*l0

def nC(x, y):
    if ((target_min < x < target_max) and (9.5*l0 < y < 25.5*l0)):
        return n_C
    else:
        return 0.

def nH(x, y):
    if ((target_min < x < target_max) and (9.5*l0 < y < 25.5*l0)):
        return n_H
    else:
        return 0.
def nE(x, y):
    if ((target_min < x < target_max) and (9.5*l0 < y < 25.5*l0)):
        return n_H+6*n_C
    else:
        return 0.

Species(
    name      = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 3,
    mass = 1.,
    atomic_number = None,
    #maximum_charge_state = None,
    number_density = nE,
    # charge_density = None,
    charge = -1.,
    mean_velocity = [0.],
    temperature = [1e-10],
    boundary_conditions = [
        ["remove", "remove"],
    #    ["periodic", "periodic"],
    #    ["periodic", "periodic"],
    ],
    # thermal_boundary_temperature = None,
    # thermal_boundary_velocity = None,
    time_frozen = 0.0,
    # ionization_model = "none",
    # ionization_electrons = None,
    # ionization_rate = None,
    is_test = False,
    pusher = "boris",

    # Radiation reaction, for particles only:
    radiation_model = "corrected-Landau-Lifshitz",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,

    # Relativistic field initialization:
    relativistic_field_initialization = "False",

    # For photon species only:
    #multiphoton_Breit_Wheeler = ["electron","positron"],
    #multiphoton_Breit_Wheeler_sampling = [1,1]

    # Merging
    merging_method = "vranic_spherical",
    merge_every = 5,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,16],
)
Species(
    name      = "proton",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 3,
    mass = 1836.,
    atomic_number = None,
    #maximum_charge_state = None,
    number_density = nH,
    # charge_density = None,
    charge = 1.,
    mean_velocity = [0.],
    temperature = [1e-10],
    boundary_conditions = [
        ["remove", "remove"],
    #    ["periodic", "periodic"],
    #    ["periodic", "periodic"],
    ],
    # thermal_boundary_temperature = None,
    # thermal_boundary_velocity = None,
    time_frozen = 0.0,
    # ionization_model = "none",
    # ionization_electrons = None,
    # ionization_rate = None,
    is_test = False,
    pusher = "boris",
# Radiation reaction, for particles only:
#     radiation_model = "corrected-Landau-Lifshitz",
#     radiation_photon_species = "photon",
#     radiation_photon_sampling = 1,
#     radiation_photon_gamma_threshold = 20,

    # Relativistic field initialization:
    relativistic_field_initialization = "False",

    # For photon species only:
    #multiphoton_Breit_Wheeler = ["electron","positron"],
    #multiphoton_Breit_Wheeler_sampling = [1,1]

    # Merging
    merging_method = "vranic_spherical",
    merge_every = 5,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,16],
)
Species(
    name      = "carbon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 3,
    mass = 1836.*12,
    atomic_number = None,
    #maximum_charge_state = None,
    number_density = nC,
    # charge_density = None,
    charge = 6.,
    mean_velocity = [0.],
    temperature = [1e-10],
    boundary_conditions = [
        ["remove", "remove"],
    #    ["periodic", "periodic"],
    #    ["periodic", "periodic"],
    ],
    # thermal_boundary_temperature = None,
    # thermal_boundary_velocity = None,
    time_frozen = 0.0,
    # ionization_model = "none",
    # ionization_electrons = None,
    # ionization_rate = None,
    is_test = False,
    pusher = "boris",
# Radiation reaction, for particles only:
#     radiation_model = "corrected-Landau-Lifshitz",
#     radiation_photon_species = "photon",
#     radiation_photon_sampling = 1,
#     radiation_photon_gamma_threshold = 20,

    # Relativistic field initialization:
    relativistic_field_initialization = "False",

    # For photon species only:
    #multiphoton_Breit_Wheeler = ["electron","positron"],
    #multiphoton_Breit_Wheeler_sampling = [1,1]

    # Merging
    merging_method = "vranic_spherical",
    merge_every = 5,
    merge_min_particles_per_cell = 16,
    merge_max_packet_size = 4,
    merge_min_packet_size = 4,
    merge_momentum_cell_size = [16,16,16],
)

RadiationReaction(
    minimum_chi_continuous = 1e-3,
#    minimum_chi_discontinuous = 1e-2,
#    table_path = "/gpfshome/mds/staff/mlobet/smilei/databases/"
)

DiagScalar(every=10*(sqrt(cpx**2+cpy**2))/0.8)
DiagFields(
    #name = "my field diag",
    every = 10*(sqrt(cpx**2+cpy**2))/0.8, #单位timestep
    #time_average = 2,
    fields = ['Ey', 'Rho_electron'],
    #subgrid = None
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 10*(sqrt(cpx**2+cpy**2))/0.8,
    time_average = 1,
    species = ["proton"],
    axes = [ ["ekin",    2,    2500,   1000, "logscale"] ]
)
