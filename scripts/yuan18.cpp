//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file yuan18.cpp
//! \brief Migrate Yuan & Yoon's work in 2018 to Athena++ platform and change to 3D.
//!
//! The formula without specific source is referenced from Ciotti & Ostriker's *AGN Feedback in Elliptical Galaxies: Numerical Simulations*

// C headers

// C++ headers
#include <sstream>  // stringstream
#include <iostream> // endl, ostream
#include <random>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../field/field.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../cooling/cooling.hpp"
#include "../eos/eos.hpp"
#include "../scalars/scalars.hpp"
#include "../units/units.hpp"
#include "../utils/buffer_utils.hpp"

// #define debug

// Debug flags
#define StellarFeedbackFlag

#ifdef StellarFeedbackFlag
#define StarEvolutionFlag
#define StarFormationFlag
#define SNiaFeedbackFlag
#ifdef StarFormationFlag
#define SNiiFeedbackFlag
#endif
#define MetalYieldsFlag
#endif
#define SN_NUM_MAX_PROC 10
#define SN_NUM_MAX_GLOBAL 200

// #define GravityFlag
#ifdef GravityFlag
// #define HaloGravityFlag
#endif
#define CheckNaN 0b00001 // 0b00001: check primitive variables, 0b00010: check conserved variables, 0b00100: check scalars, 0b01000: check user mesh block data, 0b10000: check divB constraint

#define AGNWindFlag
// #define AGNJetFlag
#if COOLING_ENABLED == 2
#define AGNRadFlag
#endif
// #define CGMInflowFlag
#define HydrostaticBCFlag

#define N_WIND 10000
#define N_JET 1000

enum UserMeshBlockIndex
{
  STAR,
  NEWSTAR,
  NOVA,       // rhoDotSNII
  NOVA_DELAY, // rhoDotSNII_DELAY
  AGN_FLUX,
  COOL_RATE,
  HEAT_RATE_AGN,
  HEAT_RATE_UVB,
  RHODOTSF,
  RHODOTEV,
  DIVB,
  N_USER_MESH_BLOCK
};
static const char *enumUserMeshBlockStr[] = {"star", "newstar", "nova", "nova_delay", "agn_flux", "cool_rate", "heat_rate_agn", "heat_rate_uvb", "rhodotsf", "rhodotev", "divB"};

// NOTE since turbulent diffusion only support working on last term of scalar
// please add a scalar before metal if needed
// to prevent confusion, please use Z for absolute metallicity and z for relative if needed
enum ScalarIndex
{
  ISM,
  CGM,
  SW,
  SNIA,
  SNII,
  AGNWC,
  AGNWH,
  AGNJ,
  METAL
};
static const char *enumScalarStr[] = {"ism_frac", "cgm_frac", "sw_frac", "snia_frac", "snii_frac", "agnwc_frac", "agnwh_frac", "agnj_frac", "metallicity"};
static const int agnhst = 7;
static const int starhst = 2;
static Real H_frac(0.7);    // dimensionless
static Real kappa_es(0.35); // in cm^2/g
static Real univAge(13.7);  // in Gyr
static const Real large_time = HUGE_NUMBER;
static Real user_timestep = HUGE_NUMBER;
static Real vmax;
static Real dt_init;
static Real redshift = 1;

// stellar feedback
static Real Mstar, rstar;
static Real d0, rc, beta;
static Real vc;
static Real etaSF;
static Real temp_sf;
static Real nrho_sf;
static Real tauDelay, tauII;
static Real rEff;
static Real sigma0;

// initial metallicity
static Real z_init_in;
static Real z_init_eff;
static Real z_init_out;
static Real z_wind;
static Real z_cgm;
static Real z_slope;
const static Real m_star[10] = {0.9, 1.0, 1.2, 1.5, 1.8, 1.9, 2.0, 2.2, 2.5, 3.0};
const static Real Z_s[6] = {0., 0.001, 0.004, 0.008, 0.02, 0.05};
const static Real yields[6][10] = {{4.57472832e-04, 1.65980683e-03, 2.56474495e-03, 2.55909551e-03,
                                    2.19045376e-03, 1.94207873e-03, 1.69500450e-03, 3.13992764e-03,
                                    4.61016491e-03, 8.01908322e-03},
                                   {0.00000000e+00, 2.76226368e-04, 3.23428480e-03, 8.38937049e-03,
                                    1.44818488e-02, 1.77502448e-02, 2.08247648e-02, 2.60002335e-02,
                                    1.88863776e-02, 8.57354671e-03},
                                   {0.00000000e+00, 6.35493925e-06, 1.18260048e-04, 1.32190537e-03,
                                    3.24840853e-03, 4.98739672e-03, 5.52185491e-03, 1.15255636e-02,
                                    1.52684984e-02, 1.41297695e-02},
                                   {0.00000000e+00, 1.21964480e-05, 2.67916024e-05, 2.18141405e-04,
                                    1.77782350e-03, 2.91128504e-03, 4.12689420e-03, 7.18149814e-03,
                                    9.16830852e-03, 1.40951338e-02},
                                   {0.00000000e+00, 2.74514687e-05, 6.20931550e-05, 8.66091012e-05,
                                    1.02388524e-04, 1.09243252e-04, 1.11612308e-04, 5.80334133e-04,
                                    2.52117768e-03, 5.89753620e-03},
                                   {0.00000000e+00, 2.74514687e-05, 6.20931550e-05, 8.66091012e-05,
                                    1.02388524e-04, 1.09243252e-04, 1.11612308e-04, 5.80334133e-04,
                                    2.52117768e-03, 5.89753620e-03}};

// gravitational potential
static Real q_star;            // axial ratio of the star profile
static Real q_galaxy;          // axial ratio of the galaxy profile
static Real eps_galaxy_radius; // r_galaxy = eps_galaxy_radius * rstar
static Real eps_galaxy_mass;   // M_galaxy = eps_galaxy_mass * Mstar
static Real v_halo;            // asymptotic circular velocity of spherically symmetric quasi-isothermal DM halo
static Real eps_halo_radius;   // r_halo = eps_halo_radius * rstar

// CGM Inflow
static Real eps_cgm_mass;
static Real t0_cgm(9.0);        // in Gyr, Ciotti 2022 ApJ
static Real timespan_cgm(12.0); // in Gyr

// MHD
static Real mhd_beta0;

// SN
static Real etaSN(0.85); // Eq 4.20
static Real massiveStar(0.1234);
static Real nrho_sn0(100); // in cm^-3
#ifdef SNiaFeedbackFlag
// Predefined SN injection radius
static Real ergSNia(1e51);    // in erg, Eq 4.11
static Real SNiaRadius0(0.9); // in pc, Martizzi 2015 mnras 450,504
static Real SNiaMass(1.4);    // in Msun
static Real hh(0.75), ss(1.1), thetaSNia(1.0);
static Real rateConst(0.22e-3); // in Gyr^(-1), This rate in Eq 4.11 is 0.32, may be a litte bigger
static Real SNiaRate0;
static Real eps_snia;
#endif // SNiaFeedbackFlag
#ifdef SNiiFeedbackFlag
// Sukhbold_2016_ApJ_821_38
// SNii progenitor ZAMS mass 9-120 Msun. During the evolution, the star loses mass.
// Only 75% of the massive star will explode, while the rest will collapse into a black hole, which will not release energy.
// We take the average of the mass change of the all massive stars as the mass of snii ejecta.
// the kinetic energy is slightly less than 1e51 erg, but we still take 1e51 erg considering the wind of newly formed massive stars and the internal energy of the ejecta.
static Real ergSNii(1e51);    // in erg
static Real mZAMS(21.34);     // in Msun, the IMF-averaged ZAMS mass of SNii progenitors
static Real SNiiRadius0(0.9); // in pc
static Real SNiiMass(16.6);   // in mSun
static Real eps_snii;
#endif // SNiiFeedbackFlag

// AGN feedback
enum UserMeshRealDataIndex
{
  M_BH,
  MDISK,
  MFALL,
  MDOT_BONDI,
  MDOT_BH,
  MDOT_WIND,
  V_WIND,
  L_BH,
  MDOT_WIND_NOW,
  PDOT_WIND_NOW,
  EPS_WIND_NOW,
  MDOT_JET_NOW,
  PDOT_JET_NOW,
  EPS_JET_NOW,
  N_USER_MESH_REAL
};
enum UserMeshIntDataIndex
{
  MODE_WIND_NOW,
  MODE_JET_NOW,
  N_USER_MESH_INT
};
static Real eff_em_factor;
static Real alpha_disk;
static Real nschr_in;
static Real nschr_out;
static Real gamma_wind = 4.0 / 3;
static Real total_area_hot;
static Real total_area_sub;
static Real total_area_sup;
static Real total_area_jet;

// for agn fb hot mode, refer to the paragraph following Yuan et al. 2018 eq.20
static Real ang1_hot;
static Real ang2_hot;
static Real ang1_sub;
static Real ang2_sub;
static Real ang1_sup;
static Real ang2_sup;
static Real ang1_jet;
static Real ang2_jet;

struct sninfo
{
  Real x1 = 0, x2 = 0, x3 = 0, radius = 0, rho = 0, vol = 0;
};
class SNInfo
{
public:
  SNInfo() {}
  void add(Real x3, Real x2, Real x1, Real radius, Real rho, Real vol)
  {
    if (num >= SN_NUM_MAX_GLOBAL)
    {
      std::stringstream msg;
      msg << "### FATAL ERROR in SNInfo::add" << std::endl
          << "Number of supernova exceeds the maximum number " << SN_NUM_MAX_GLOBAL << std::endl
          << "Please increase SN_NUM_MAX_GLOBAL if you think it is right to have so many SN in ONE timestep" << std::endl;
      ATHENA_ERROR(msg);
    }
    info[num].x1 = x1;
    info[num].x2 = x2;
    info[num].x3 = x3;
    info[num].radius = radius;
    info[num].rho = rho;
    info[num].vol = vol;
    num++;
  }
  void clear()
  {
    num = 0;
  }
  sninfo info[SN_NUM_MAX_GLOBAL];
  int num = 0;
} SNiaInfo, SNiiInfo;

enum OutflowMode
{
  NONE,
  HOT,
  SUB,
  SUP,
  JET
};
enum OutflowInfo
{
  MDOT,
  PDOT,
  EPS_TE,
  TIMESPAN,
  TIME_ARRIVAL,
};
class Outflow
{
public:
  Outflow(AthenaArray<Real> &info0, AthenaArray<int> &mode0)
  {
    info = &info0;
    mode = &mode0;
    n_ = mode0.GetSize();
  }

  void pop_front(Real time_now, Real timestep, Real &mdot_now, Real &pdot_now, Real &eps_te_now, int &mode_now)
  {
    mdot_now = 0;
    pdot_now = 0;
    eps_te_now = 0;
    mode_now = NONE;
    Real timespan = 0;
    Real timespan_left = 0;
    // accumulate the outflow that will pass the inner boundary
    int i;
    for (i = 0; i < n_; i++)
    {
      if (time_now < (*info)(i, TIME_ARRIVAL))
      {
        break;
      }
      mdot_now += (*info)(i, MDOT) * (*info)(i, TIMESPAN);
      pdot_now += (*info)(i, PDOT) * (*info)(i, TIMESPAN);
      eps_te_now = (*info)(i, EPS_TE);
      timespan = (*info)(i, TIMESPAN);
      mode_now = (*mode)(i);
    }
    if (timespan > TINY_NUMBER)
    {
      mdot_now /= timespan;
      pdot_now /= timespan;

      bool is_left = false;
      if (timespan > timestep)
      {
        timespan_left = timespan - timestep;
        (*info)(0, MDOT) = mdot_now;
        (*info)(0, PDOT) = pdot_now;
        (*info)(0, EPS_TE) = eps_te_now;
        (*info)(0, TIMESPAN) = timespan_left;
        (*info)(0, TIME_ARRIVAL) = time_now + timestep;
        (*mode)(0) = mode_now;
        is_left = true;
      }
      else
      {
        mdot_now *= timespan / timestep;
        pdot_now *= timespan / timestep;
      }
      // remove the outflow that has passed the inner boundary
      int j = is_left;
      for (; j < n_ - i; j++)
      {
        if ((*mode)(j + i - is_left) == NONE)
        {
          break;
        }
        (*info)(j, MDOT) = (*info)(j + i - is_left, MDOT);
        (*info)(j, PDOT) = (*info)(j + i - is_left, PDOT);
        (*info)(j, EPS_TE) = (*info)(j + i - is_left, EPS_TE);
        (*info)(j, TIMESPAN) = (*info)(j + i - is_left, TIMESPAN);
        (*info)(j, TIME_ARRIVAL) = (*info)(j + i - is_left, TIME_ARRIVAL);
        (*mode)(j) = (*mode)(j + i - is_left);
      }
      // clear the rest of the buffer
      for (int k = j; k < j + i - is_left; k++)
      {
        (*info)(k, MDOT) = 0;
        (*info)(k, PDOT) = 0;
        (*info)(k, EPS_TE) = 0;
        (*info)(k, TIMESPAN) = 0;
        (*info)(k, TIME_ARRIVAL) = large_time;
        (*mode)(k) = NONE;
      }
    }
  };
  void push_back(Real mdot, Real pdot, Real eps_te, Real timespan, Real time_arrival, int mode_flag)
  {
    // find the position to insert
    int i;
    Real v_travel = pdot / mdot;
    if (time_arrival > univAge)
      time_arrival = univAge;
    for (i = 0; i < n_; i++)
    {
      if (time_arrival <= (*info)(i, TIME_ARRIVAL))
      {
        break;
      }
    }
    if (i == n_)
    {
      std::stringstream msg;
      msg << "### FATAL ERROR in Outflow::push_back" << std::endl
          << "Outflow buffer is full, please increase N_WIND or N_JET" << std::endl;
      ATHENA_ERROR(msg);
    }
    // accumulate the wind that will be pushed by the new wind and remove them
    for (int j = i; j < n_; j++)
    {
      if ((*mode)(j) == NONE)
      {
        break;
      }
      mdot += (*info)(j, MDOT) * (*info)(j, TIMESPAN) / timespan;
      pdot += (*info)(j, PDOT) * (*info)(j, TIMESPAN) / timespan;
      (*info)(j, MDOT) = 0;
      (*info)(j, PDOT) = 0;
      (*info)(j, EPS_TE) = 0;
      (*info)(j, TIMESPAN) = 0;
      (*info)(j, TIME_ARRIVAL) = large_time;
      (*mode)(j) = NONE;
    }
    // insert the new wind
    (*info)(i, MDOT) = mdot;
    (*info)(i, PDOT) = pdot;
    (*info)(i, EPS_TE) = eps_te;
    (*info)(i, TIMESPAN) = timespan;
    (*info)(i, TIME_ARRIVAL) = time_arrival;
    (*mode)(i) = mode_flag;
  };

private:
  AthenaArray<Real> *info;
  AthenaArray<int> *mode;
  int n_;
};

//
// SECTION Function Declarations
void yuan18Source(MeshBlock *, const Real, const Real, const AthenaArray<Real> &, const AthenaArray<Real> &, const AthenaArray<Real> &, AthenaArray<Real> &, AthenaArray<Real> &);
void InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
Real UserTimeStep(MeshBlock *pmb);
Real MeshGeneratorX2(Real x, RegionSize rs);

Real expFactor(Real dt, Real tau);

Real GetTotalAcceleration(Real r);
Real PotentialFromGalaxy(Real r);
Real AccelerationFromGalaxy(Real r);
#ifdef GravityFlag
Real PotentialFromBH(Real r);
Real AccelerationFromBH(Real r);
#ifdef HaloGravityFlag
Real PotentialFromHalo(Real r);
Real AccelerationFromHalo(Real r);
#endif // HaloGravityFlag
#endif // GravityFlag

#ifdef CheckNaN
bool CheckCondition(PrimIndex, Real);
bool CheckCondition(ConsIndex, Real);
bool CheckCondition(ScalarIndex, Real);
bool CheckCondition(UserMeshBlockIndex, Real);
#endif

Real GetMdot_Bondi(Mesh *pm, Coordinates *pco, const AthenaArray<Real> &prim, int is, int js, int je, int ks, int ke);
Real GetMdot_BH_Cold(Real mdot_in);
Real MdotWind_Cold(Real mdot_bh);
int Newtonian_Solver(Real x0, Real eps, const std::function<Real(Real)> &f, Real *root);
Real AGNHistoryOutput(MeshBlock *pmb, int iout);
Real StarHistoryOutput(MeshBlock *pmb, int iout);

Real GetStarEvoMetal(Real mTO_sun, Real z_star);
Real GetStarMetal(Real r);
Real GetRadEfficiency(Real mdot);
Real GetMeanMolWeight(Real Z_Zsun);
Real GetTemp(Real press, Real rho, Real Z_Zsun);
Real GetPress(Real temp, Real rho, Real Z_Zsun);

// #SECTION Function Declarations
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Enroll AGN boundary and sourceterm
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  kappa_es *= SQR(punit->cm_code) / punit->gram_code;
  univAge = 13.7 * punit->giga_yr_code;

  dt_init = pin->GetReal("problem", "dt_init") * punit->yr_code;
  vmax = pin->GetReal("problem", "vmax") * punit->speed_of_light_code;

  Mstar = pin->GetReal("problem", "Mstar") * punit->solar_mass_code;
  d0 = pin->GetReal("problem", "d0") * punit->hydrogen_mass_code / CUBE(punit->cm_code);
  rc = pin->GetReal("problem", "rc") * punit->kpc_code;
  // beta = pin->GetReal("problem", "beta");// TODO If it changes, it needs to be re-integrated
  beta = 2. / 3.;

  etaSF = pin->GetReal("problem", "etaSF");
  temp_sf = pin->GetOrAddReal("problem", "temp_sf", 4e4) * punit->kelvin_code;
  nrho_sf = pin->GetOrAddReal("problem", "nrho_sf", 1) / CUBE(punit->cm_code);
  tauDelay = pin->GetReal("problem", "tauDelay") * punit->giga_yr_code;
  tauII = pin->GetReal("problem", "tauII") * punit->giga_yr_code;
  sigma0 = pin->GetReal("problem", "sigma0") * punit->km_s_code;

  z_init_in = pin->GetOrAddReal("problem", "z_init_in", 2);
  z_init_eff = pin->GetOrAddReal("problem", "z_init_eff", 0.9);
  z_cgm = pin->GetOrAddReal("problem", "z_cgm", 0.1);
  z_wind = pin->GetOrAddReal("problem", "z_wind", 1.4);
  z_slope = pin->GetOrAddReal("problem", "z_slope", -0.23);

#ifndef GravityFlag
  rstar = pin->GetReal("problem", "rstar") * punit->kpc_code;
  rEff = pin->GetReal("problem", "rEff") * punit->kpc_code;
  vc = pin->GetReal("problem", "vc") * punit->km_s_code;
#else
  eps_galaxy_radius = pin->GetReal("problem", "eps_galaxy_radius");
  q_star = pin->GetReal("problem", "q_star");
  q_galaxy = pin->GetReal("problem", "q_galaxy");
  eps_galaxy_mass = eps_galaxy_radius / q_star;
  rstar = rc / (0.7447 * std::sqrt(q_star));
#endif
#ifdef HaloGravityFlag
  v_halo = pin->GetReal("problem", "v_halo") * punit->km_s_code;
  eps_halo_radius = pin->GetReal("problem", "eps_halo_radius");
#endif

#ifdef CGMInflowFlag
  eps_cgm_mass = pin->GetReal("problem", "eps_cgm_mass");
  t0_cgm *= punit->Gyr_code;
  timespan_cgm *= punit->Gyr_code;
#endif

  ang1_hot = pin->GetOrAddReal("problem", "ang1_hot", 30.0) / 180.0 * PI;
  ang2_hot = pin->GetOrAddReal("problem", "ang2_hot", 70.0) / 180.0 * PI;
  ang1_sub = pin->GetOrAddReal("problem", "ang1_sub", 0.0) / 180.0 * PI;
  ang2_sub = pin->GetOrAddReal("problem", "ang2_sub", 0.0) / 180.0 * PI;
  ang1_sup = pin->GetOrAddReal("problem", "ang1_sup", 0.0) / 180.0 * PI;
  ang2_sup = pin->GetOrAddReal("problem", "ang2_sup", 30.0) / 180.0 * PI;
  ang1_jet = pin->GetOrAddReal("problem", "ang1_jet", 0.0) / 180.0 * PI;
  ang2_jet = pin->GetOrAddReal("problem", "ang2_jet", 10.0) / 180.0 * PI;
  // MHD
  if (MAGNETIC_FIELDS_ENABLED)
  {
    mhd_beta0 = pin->GetOrAddReal("problem", "mhd_beta0", 10.0);
  }

  nrho_sn0 *= punit->hydrogen_mass_code / CUBE(punit->cm_code);
#ifdef SNiaFeedbackFlag
  ergSNia *= punit->erg_code;
  SNiaRadius0 *= punit->pc_code;
  SNiaMass *= punit->solar_mass_code;
  rateConst /= punit->giga_yr_code;
  SNiaRate0 = rateConst * SQR(hh) * thetaSNia;
  eps_snia = etaSN * ergSNia / SNiaMass;
#endif // SNiaFeedbackFlag
#ifdef SNiiFeedbackFlag
  ergSNii *= punit->erg_code;
  mZAMS *= punit->solar_mass_code;
  SNiiRadius0 *= punit->pc_code;
  SNiiMass *= punit->solar_mass_code;
  eps_snii = etaSN * ergSNii / SNiiMass;
#endif // SNiiFeedbackFlag
  const int nhst = agnhst + starhst;
  AllocateUserHistoryOutput(nhst);
  const char *AGNHistoryVariableName[agnhst] = {"m_bh", "mdot_in", "mdot_bh", "mdot_wind", "v_wind", "l_bh", "mode_now"};
  const char *StarHistoryVariableName[starhst] = {"star", "sfr"};
  for (int i = 0; i < agnhst; i++)
  {
    EnrollUserHistoryOutput(i, AGNHistoryOutput, AGNHistoryVariableName[i], UserHistoryOperation::max);
  }
  for (int i = 0; i < starhst; i++)
  {
    EnrollUserHistoryOutput(agnhst + i, StarHistoryOutput, StarHistoryVariableName[i], UserHistoryOperation::sum);
  }

  AllocateRealUserMeshDataField(3);
  ruser_mesh_data[0].NewAthenaArray(N_USER_MESH_REAL);
  ruser_mesh_data[1].NewAthenaArray(N_WIND, 5);
  ruser_mesh_data[2].NewAthenaArray(N_JET, 5);
  AllocateIntUserMeshDataField(3);
  iuser_mesh_data[0].NewAthenaArray(N_USER_MESH_INT);
  iuser_mesh_data[1].NewAthenaArray(N_WIND);
  iuser_mesh_data[2].NewAthenaArray(N_JET);
  for (int i = 0; i < N_WIND; i++)
  {
    ruser_mesh_data[1](i, MDOT) = 0;
    ruser_mesh_data[1](i, PDOT) = 0;
    ruser_mesh_data[1](i, EPS_TE) = 0;
    ruser_mesh_data[1](i, TIMESPAN) = 0;
    ruser_mesh_data[1](i, TIME_ARRIVAL) = large_time;
    iuser_mesh_data[1](i) = NONE;
  }
  for (int i = 0; i < N_JET; i++)
  {
    ruser_mesh_data[2](i, MDOT) = 0;
    ruser_mesh_data[2](i, PDOT) = 0;
    ruser_mesh_data[2](i, EPS_TE) = 0;
    ruser_mesh_data[2](i, TIMESPAN) = 0;
    ruser_mesh_data[2](i, TIME_ARRIVAL) = large_time;
    iuser_mesh_data[2](i) = NONE;
  }
  ruser_mesh_data[0](M_BH) = pin->GetReal("problem", "eps_bh_mass") * Mstar;
  eff_em_factor = pin->GetReal("problem", "eff_em_factor");
  alpha_disk = pin->GetReal("problem", "alpha_disk");
  nschr_in = pin->GetReal("problem", "nschr_in");

  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user"))
  {
    EnrollUserBoundaryFunction(inner_x1, InnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user"))
  {
    EnrollUserBoundaryFunction(outer_x1, OuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user"))
  {
    EnrollUserBoundaryFunction(inner_x2, InnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user"))
  {
    EnrollUserBoundaryFunction(outer_x2, OuterX2);
  }
  EnrollUserExplicitSourceFunction(yuan18Source); // TODO may need move some sourceterm to src/hydro/srcterms
  EnrollUserTimeStepFunction(UserTimeStep);
  if (pin->GetReal("mesh", "x2rat") < 0.0)
  {
    EnrollUserMeshGenerator(X2DIR, MeshGeneratorX2);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes gas and star density.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
#ifdef AGNRadFlag
  // first, ensure that # of processors = # of mesh blocks which is required by the AGN feedback module
  if (Globals::nranks != pmy_mesh->nbtotal)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "Number of processors must be equal to the number of mesh blocks "
        << "(required by the AGN feedback module), i.e.," << std::endl
        << "(mesh_size.nx1 / block_size.nx1) * (mesh_size.nx2 / block_size.nx2) "
           "* (mesh_size.nx3 / block_size.nx3) must be equal to the number of processors"
        << std::endl
        << "### Current configuration:" << std::endl
        << "mesh_size.nx1 = " << pmy_mesh->mesh_size.nx1 << ", mesh_size.nx2 = " << pmy_mesh->mesh_size.nx2
        << ", mesh_size.nx3 = " << pmy_mesh->mesh_size.nx3 << std::endl
        << "block_size.nx1 = " << block_size.nx1 << ", block_size.nx2 = " << block_size.nx2 << ", block_size.nx3 = "
        << block_size.nx3 << std::endl
        << "Number of mesh blocks: " << pmy_mesh->nbtotal << std::endl
        << "Number of processors: " << Globals::nranks << std::endl;
    ATHENA_ERROR(msg);
  }
#endif

  pmy_mesh->dt = dt_init;
  Real x1, x2, x3;
  Real starDen(0.0), gasDen(0.0);
  AthenaArray<Real> press(ncells3, ncells2, ncells1);
  Real gammaGas = peos->GetGamma();
  Real gamma1 = gammaGas - 1.0;

  for (int k = ks; k <= ke; ++k)
  {
    for (int j = js; j <= je; ++j)
    {
      for (int i = is; i <= ie; ++i)
      {
        x1 = pcoord->x1v(i);
        for (int n = 0; n < NSCALARS; ++n)
        {
          pscalars->s(n, k, j, i) = 0.0;
        }
        // Gas setup
        Real ratio1 = (rc * rc) / (x1 * x1);
#ifndef GravityFlag

        press(k, j, i) = 0.5 * d0 * vc * vc * (std::log(1 + ratio1) - ratio1 / (1 + ratio1)); // Integrate[rho[x]*g[x],{x,r,Infinity}]

#else

        Real rgal = eps_galaxy_radius * rstar;
        Real ratio2 = SQR(rgal / rc);
        Real ratio3 = rgal / x1;
        press(k, j, i) = d0 * G * Mstar * eps_galaxy_mass *
                             (std::log(1 + ratio1) / (2 * rgal) +
                              std::log(SQR(ratio3 + 1) / (ratio1 + 1)) / (2 * rgal * SQR(1 + ratio2)) +
                              (3 + ratio2) / (4 * rc * SQR(1 + ratio2)) * (2 * std::atan(x1 / rc) - PI) +
                              (x1 - rgal) / (2 * (SQR(rc) + SQR(x1)) * (1 + ratio2)))                         // Galaxy
                         + d0 * G * m_bh * (2 / x1 + x1 / (SQR(x1) + SQR(rc)) - 3 * std::atan(rc / x1) / rc); // Black hole

#ifdef HaloGravityFlag

        // TODO update this press by adding halo term.

#endif                                                     // HaloGravityFlag
#endif                                                     // GravityFlag
        gasDen = d0 * std::pow(1 + 1 / ratio1, -3 * beta); // Gas density given by deprojected beta profile

        phydro->u(IDN, k, j, i) = gasDen;
        phydro->u(IM1, k, j, i) = 0.0;
        phydro->u(IM2, k, j, i) = 0.0;
        phydro->u(IM3, k, j, i) = 0.0;
        phydro->u(IEN, k, j, i) = press(k, j, i) / gamma1;
        phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM1, k, j, i)) / phydro->u(IDN, k, j, i);
        phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM2, k, j, i)) / phydro->u(IDN, k, j, i);
        phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM3, k, j, i)) / phydro->u(IDN, k, j, i);

        phydro->w(IDN, k, j, i) = gasDen;
        phydro->w(IVX, k, j, i) = phydro->u(IM1, k, j, i) / phydro->u(IDN, k, j, i);
        phydro->w(IVY, k, j, i) = phydro->u(IM2, k, j, i) / phydro->u(IDN, k, j, i);
        phydro->w(IVZ, k, j, i) = phydro->u(IM3, k, j, i) / phydro->u(IDN, k, j, i);
        phydro->w(IPR, k, j, i) = press(k, j, i);

        if (NSCALARS > 0)
        {
          pscalars->s(ISM, k, j, i) = gasDen; // Conservative form of scalars
          // Metallicity setup
          pscalars->s(METAL, k, j, i) = gasDen * GetStarMetal(x1) * Constants::solar_metallicity;
        }
// End gas setup

// Star setup
#ifndef GravityFlag
        starDen = (Mstar * rstar) / (4 * PI * SQR(x1 * (x1 + rstar))); // Yoon et al. ApJ 864:6, 2018, Eq(3)
#else
        Real hstar = x1 / rstar * std::sqrt(SQR(std::cos(x2)) + SQR(std::sin(x2) / q_star));
        starDen = Mstar / (4 * PI * q_star * CUBE(rstar) * SQR(hstar * (1 + hstar)));
#endif

        ruser_meshblock_data[STAR](k, j, i) = starDen;
        // End star setup
      }
    }
  }

  // Magnetic field setup
  if (MAGNETIC_FIELDS_ENABLED)
  {
    int nx1 = ie - is + 1 + 2 * NGHOST;
    int nx2 = 1;
    if (je > js)
      nx2 = je - js + 1 + 2 * NGHOST;
    int nx3 = 1;
    if (ke > ks)
      nx3 = ke - ks + 1 + 2 * NGHOST;
    EdgeField mvp(nx3, nx2, nx1); // magnetic vector potential
    AthenaArray<Real> area(nx1);
    for (int k = ks; k <= ke + 1; ++k)
    {
      for (int j = js; j <= je + 1; ++j)
      {
        for (int i = is; i <= ie; ++i)
        {
          mvp.x1e(k, j, i) = 0;
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        for (int i = is; i <= ie + 1; ++i)
        {
          mvp.x2e(k, j, i) = 0;
        }
      }
    }
    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je + 1; ++j)
      {
        for (int i = is; i <= ie + 1; ++i)
        {
          mvp.x3e(k, j, i) = std::sqrt(2 * press(k, j, i) / (mhd_beta0 * std::sqrt(pcoord->x1f(i) / (1.0 * punit->pc_code)))) * pcoord->x1f(i) * std::sin(pcoord->x2f(j)); // magnetic field magnitude will be sqrt(2P/beta) at z axis;
        }
      }
    }
    if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
    {
      for (int k = ks; k <= ke; ++k)
      {
        for (int j = js; j <= je; ++j)
        {
          pcoord->Face1Area(k, j, is, ie + 1, area);
          for (int i = is; i <= ie + 1; ++i)
          {
            pfield->b.x1f(k, j, i) = ((pcoord->GetEdge3Length(k, j + 1, i) * mvp.x3e(k, j + 1, i) - pcoord->GetEdge3Length(k, j, i) * mvp.x3e(k, j, i)) - (pcoord->GetEdge2Length(k + 1, j, i) * mvp.x2e(k + 1, j, i) - pcoord->GetEdge2Length(k, j, i) * mvp.x2e(k, j, i))) / area(i);
          }
        }
      }
      for (int k = ks; k <= ke; ++k)
      {
        for (int j = js; j <= je + 1; ++j)
        {
          pcoord->Face2Area(k, j, is, ie, area);
          for (int i = is; i <= ie; ++i)
          {
            pfield->b.x2f(k, j, i) = ((pcoord->GetEdge1Length(k + 1, j, i) * mvp.x1e(k + 1, j, i) - pcoord->GetEdge1Length(k, j, i) * mvp.x1e(k, j, i)) - (pcoord->GetEdge3Length(k, j, i + 1) * mvp.x3e(k, j, i + 1) - pcoord->GetEdge3Length(k, j, i) * mvp.x3e(k, j, i))) / area(i);
          }
        }
      }
      for (int k = ks; k <= ke + 1; ++k)
      {
        for (int j = js; j <= je; ++j)
        {
          pcoord->Face3Area(k, j, is, ie, area);
          for (int i = is; i <= ie; ++i)
          {
            pfield->b.x3f(k, j, i) = ((pcoord->GetEdge2Length(k, j, i + 1) * mvp.x2e(k, j, i + 1) - pcoord->GetEdge2Length(k, j, i) * mvp.x2e(k, j, i)) - (pcoord->GetEdge1Length(k, j + 1, i) * mvp.x1e(k, j + 1, i) - pcoord->GetEdge1Length(k, j, i) * mvp.x1e(k, j, i))) / area(i);
            // pfield->b.x3f(k, j, i) += std::sqrt(2 * press(k, j, i) / (mhd_beta0 * std::sqrt(pcoord->x1f(i) / (1.0 * punit->pc_code))));
          }
        }
      }
      pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks, ke);
      for (int k = ks; k <= ke; ++k)
      {
        for (int j = js; j <= je; ++j)
        {
          for (int i = is; i <= ie; ++i)
          {
            phydro->u(IEN, k, j, i) += 0.5 * (SQR(pfield->bcc(IB1, k, j, i)) +
                                              SQR(pfield->bcc(IB2, k, j, i)) +
                                              SQR(pfield->bcc(IB3, k, j, i)));
          }
        }
      }
      for (int k = ks; k <= ke; k++)
      {
        for (int j = js; j <= je; j++)
        {
          for (int i = is; i <= ie; i++)
          {
            Real div = 0, b = 0, ratio = 0, dr = 0;
            div += pcoord->GetFace1Area(k, j, i + 1) * pfield->b.x1f(k, j, i + 1) - pcoord->GetFace1Area(k, j, i) * pfield->b.x1f(k, j, i);
            div += pcoord->GetFace2Area(k, j + 1, i) * pfield->b.x2f(k, j + 1, i) - pcoord->GetFace2Area(k, j, i) * pfield->b.x2f(k, j, i);
            div += pcoord->GetFace3Area(k + 1, j, i) * pfield->b.x3f(k + 1, j, i) - pcoord->GetFace3Area(k, j, i) * pfield->b.x3f(k, j, i);
            div /= pcoord->GetCellVolume(k, j, i);
            b = std::sqrt(SQR(pfield->bcc(IB1, k, j, i)) + SQR(pfield->bcc(IB2, k, j, i)) + SQR(pfield->bcc(IB3, k, j, i)));
            dr = std::sqrt(3. / (1 / SQR(pcoord->dx1f(i)) + 1 / SQR(pcoord->x1v(i) * pcoord->dx2f(j)) + 1 / SQR(pcoord->x1v(i) * std::sin(pcoord->x2v(j)) * pcoord->dx3f(k))));
            ratio = std::abs(div) / (b / dr);
            ruser_meshblock_data[DIVB](k, j, i) = ratio;
            if (ratio > 1e-14)
            {
              std::stringstream msg;
              msg << "\n### B FIELD ERROR ###\n"
                  << "  in cell (k,j,i) = (" << k << "," << j << "," << i << ")\n"
                  << "  found div B = " << div << ", B/dr = " << (b / dr) << ", ratio = " << ratio << "\n"
                  << "  div B should be zero!\n"
                  << "### PROGRAM ABORT ###\n";
              ATHENA_ERROR(msg);
            }
          }
        }
      }
    }
    else
    {
      std::stringstream msg;
      msg << "### FATAL ERROR in Problem Generator ###" << std::endl
          << "Magnetic field setup only works for spherical_polar coordinate system" << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}

// MeshBlock::InitUserMeshBlockData() is called before MeshBlock::ProblemGenerator()
// while MeshBlock::ProblemGenerator() is not called when restart.
// It is safe to initialize user local data which is not necessary to be stored in UserMeshData or UserMeshBlockData in InitUserMeshBlockData().
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(N_USER_MESH_BLOCK);
  AllocateRealUserMeshBlockDataField(N_USER_MESH_BLOCK);
  for (int n = 0; n < N_USER_MESH_BLOCK; n++)
  {
    SetUserOutputVariableName(n, enumUserMeshBlockStr[n]);
    ruser_meshblock_data[n].NewAthenaArray(ncells3, ncells2, ncells1);
  }
  total_area_hot = 0;
  total_area_sub = 0;
  total_area_sup = 0;
  total_area_jet = 0;
  if (pcoord->x1v(is) < pmy_mesh->mesh_size.x1min + pcoord->dx1f(is))
  {
    for (int k = ks; k <= ke; k++)
    {
      for (int j = js; j <= je; j++)
      {
        Real theta = pcoord->x2v(j);
        Real area = pcoord->GetFace1Area(k, j, is);
        if ((theta > ang1_hot && theta < ang2_hot) || (theta > PI - ang2_hot && theta < PI - ang1_hot))
        {
          total_area_hot += area;
        }
        if (ang2_sub > TINY_NUMBER)
        {
          if ((theta > ang1_sub && theta < ang2_sub) || (theta > PI - ang2_sub && theta < PI - ang1_sub))
          {
            total_area_sub += area;
          }
        }
        if ((theta > ang1_sup && theta < ang2_sup) || (theta > PI - ang2_sup && theta < PI - ang1_sup))
        {
          total_area_sup += area;
        }
        if ((theta > ang1_jet && theta < ang2_jet) || (theta > PI - ang2_jet && theta < PI - ang1_jet))
        {
          total_area_jet += area;
        }
      }
    }
  }
#ifdef MPI_PARALLEL
  Real total_area_hot_global = 0;
  Real total_area_sub_global = 0;
  Real total_area_sup_global = 0;
  Real total_area_jet_global = 0;
  MPI_Allreduce(&total_area_hot, &total_area_hot_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_area_sub, &total_area_sub_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_area_sup, &total_area_sup_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_area_jet, &total_area_jet_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  total_area_hot = total_area_hot_global;
  total_area_sub = total_area_sub_global;
  total_area_sup = total_area_sup_global;
  total_area_jet = total_area_jet_global;
#endif // MPI_PARALLEL
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for (int n = 0; n < N_USER_MESH_BLOCK; n++)
  {
    for (int k = ks; k <= ke; k++)
    {
      for (int j = js; j <= je; j++)
      {
        for (int i = is; i <= ie; i++)
        {
          user_out_var(n, k, j, i) = ruser_meshblock_data[n](k, j, i);
        }
      }
    }
  }
  for (int n = 0; n < NSCALARS; n++)
  {
    user_scalar_names_prim_[n] = enumScalarStr[n];
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkBeforeCycle()
//! \brief Function called in main before each cycle.
//========================================================================================
void Mesh::UserWorkBeforeCycle()
{
  user_timestep = std::numeric_limits<Real>::max();
  for (int nb = 0; nb < nblocal; nb++)
  {
    MeshBlock *pmb = my_blocks(nb);
    Coordinates *pcoord = pmb->pcoord;
    int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
    const AthenaArray<Real> &prim = pmb->phydro->w;
    const AthenaArray<Real> &psr = pmb->pscalars->r;
#ifdef AGNWindFlag
    Real &mdot_wind_now = pmb->pmy_mesh->ruser_mesh_data[0](MDOT_WIND_NOW);
    Real &pdot_wind_now = pmb->pmy_mesh->ruser_mesh_data[0](PDOT_WIND_NOW);
    Real &eps_wind_now = pmb->pmy_mesh->ruser_mesh_data[0](EPS_WIND_NOW);
    int &mode_wind_now = pmb->pmy_mesh->iuser_mesh_data[0](MODE_WIND_NOW);
    Outflow wind(pmb->pmy_mesh->ruser_mesh_data[1], pmb->pmy_mesh->iuser_mesh_data[1]);
    // calculate the wind coming to the inner boundary in next loop
    wind.pop_front(time, dt, mdot_wind_now, pdot_wind_now, eps_wind_now, mode_wind_now);
#endif
#ifdef AGNJetFlag
    Real &mdot_jet_now = pmb->pmy_mesh->ruser_mesh_data[0](MDOT_JET_NOW);
    Real &pdot_jet_now = pmb->pmy_mesh->ruser_mesh_data[0](PDOT_JET_NOW);
    Real &eps_jet_now = pmb->pmy_mesh->ruser_mesh_data[0](EPS_JET_NOW);
    int &mode_jet_now = pmb->pmy_mesh->iuser_mesh_data[0](MODE_JET_NOW);
    Outflow jet(pmb->pmy_mesh->ruser_mesh_data[2], pmb->pmy_mesh->iuser_mesh_data[2]);
    // calculate the jet coming to the inner boundary in next loop
    jet.pop_front(time, dt, mdot_jet_now, pdot_jet_now, eps_jet_now, mode_jet_now);
#endif // AGNJetFlag
#if defined SNiaFeedbackFlag || defined SNiiFeedbackFlag
    // std::random_device seed; // This will generate a random seed at each meshblock because each process has its own meshblock
    long seed = nbtotal * ncycle + pmb->gid;
    std::mt19937_64 gen(seed);
#ifdef SNiaFeedbackFlag
    AthenaArray<Real> &star = pmb->ruser_meshblock_data[STAR];
    // SNiaProp((k_coord,j_coord,i_corrd,snradius,sn_rho,total_vol), proc*num)
    AthenaArray<Real> SNiaProp_l(Globals::nranks * SN_NUM_MAX_PROC, 6), SNiaProp_g(Globals::nranks * SN_NUM_MAX_PROC, 6);
    int SNia_num_perproc = 0;
#endif // SNiaFeedbackFlag
#ifdef SNiiFeedbackFlag
    AthenaArray<Real> &nova = pmb->ruser_meshblock_data[NOVA];
    AthenaArray<Real> SNiiProp_l(Globals::nranks * SN_NUM_MAX_PROC, 6), SNiiProp_g(Globals::nranks * SN_NUM_MAX_PROC, 6);
    int SNii_num_perproc = 0;
#endif // SNiiFeedbackFlag

    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        for (int i = is; i <= ie; ++i)
        {

          Real x3 = pcoord->x3v(k);
          Real x2 = pcoord->x2v(j);
          Real x1 = pcoord->x1v(i);
          Real vol = pcoord->GetCellVolume(k, j, i);
          Real zmetal = psr(METAL, k, j, i) / Constants::solar_metallicity;

#ifdef SNiaFeedbackFlag
          Real starMass = star(k, j, i) * vol;
          Real Lb = starMass / (5.81 * punit->solar_mass_code);                // In Lsun
          Real SNiaDotExpect = SNiaRate0 * Lb * std::pow(time / univAge, -ss); // Eq 4.11
          std::poisson_distribution<int> SNiaPoissonDist(SNiaDotExpect * dt);
          int SNiaNum = SNiaPoissonDist(gen);
          Real SNiaRadius = SNiaRadius0 * std::pow(prim(IDN, k, j, i) * H_frac / nrho_sn0, -0.33) * std::pow(zmetal, 0.046);
          if (SNia_num_perproc >= SN_NUM_MAX_PROC)
          {
            printf("Too many SNIa!\n");
          }
          else
          {
            if (SNiaNum > 0)
            {
              SNiaProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNia_num_perproc, 0) = x3;
              SNiaProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNia_num_perproc, 1) = x2;
              SNiaProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNia_num_perproc, 2) = x1;
              SNiaProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNia_num_perproc, 3) = SNiaRadius;
              SNiaProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNia_num_perproc, 4) = SNiaNum * SNiaMass;
              SNia_num_perproc += 1;
              star(k, j, i) -= SNiaNum * SNiaMass / vol;
#ifdef debug
              printf("SNiaNum: %d at %e kpc, %e deg, %e deg, r = %e pc, mass = %e Msun\n", SNiaNum, x1 / punit->kpc_code, x2 / PI * 180, x3 / PI * 180, SNiaRadius / punit->pc_code, SNiaNum * SNiaMass / punit->solar_mass_code);
#endif
            }
          }
#endif // SNiaFeedbackFlag

#ifdef SNiiFeedbackFlag
          Real SNiiExpect = nova(k, j, i) * vol * dt / mZAMS; // nova is used to store the ZAMS mass of the star which will become snii.
          std::poisson_distribution<int> SNiiPoissonDist(SNiiExpect);
          int SNiiNum = SNiiPoissonDist(gen);
          Real SNiiRadius = SNiiRadius0 * std::pow(prim(IDN, k, j, i) * H_frac / nrho_sn0, -0.33) * std::pow(zmetal, 0.046);
          if (SNii_num_perproc >= SN_NUM_MAX_PROC)
          {
            printf("Too many SNII!\n");
          }
          else
          {
            if (SNiiNum > 0)
            {
              SNiiProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNii_num_perproc, 0) = x3;
              SNiiProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNii_num_perproc, 1) = x2;
              SNiiProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNii_num_perproc, 2) = x1;
              SNiiProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNii_num_perproc, 3) = SNiiRadius;
              SNiiProp_l(Globals::my_rank * SN_NUM_MAX_PROC + SNii_num_perproc, 4) = SNiiNum * SNiiMass;
              SNii_num_perproc += 1;
              nova(k, j, i) -= SNiiNum * mZAMS / vol; // The remnant mass should be subtracted, while ejecta mass will be added to gas mass.
#ifdef debug
              printf("SNiiNum: %d at %e kpc, %e deg, %e deg\n", SNiiNum, x1, x2 / PI * 180, x3 / PI * 180);
#endif
            }
          }
#endif // SNiiFeedbackFlag
        }
      }
    }
#ifdef SNiaFeedbackFlag
    // MPI broadcast SNia location
    MPI_Allreduce(SNiaProp_l.data(), SNiaProp_g.data(), 6 * Globals::nranks * SN_NUM_MAX_PROC, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    // Compute SNia volume
    for (int n = 0; n < Globals::nranks * SN_NUM_MAX_PROC; ++n)
    {
      Real sn_x1, sn_x2, sn_x3, sn_radius;
      sn_x3 = SNiaProp_g(n, 0);
      sn_x2 = SNiaProp_g(n, 1);
      sn_x1 = SNiaProp_g(n, 2);
      sn_radius = SNiaProp_g(n, 3);
      if (sn_radius > TINY_NUMBER)
      {
        for (int k = ks; k <= ke; ++k)
        {
          for (int j = js; j <= je; ++j)
          {
            for (int i = is; i <= ie; ++i)
            {
              Real x1 = pcoord->x1v(i), x2 = pcoord->x2v(j), x3 = pcoord->x3v(k);
              Real distance2 = SQR(x1) + SQR(sn_x1) - 2 * x1 * sn_x1 * (std::cos(x2) * std::cos(sn_x2) + std::cos(x3 - sn_x3) * std::sin(x2) * std::sin(sn_x2));
              if (distance2 < SQR(sn_radius))
              {
                SNiaProp_l(n, 5) += pcoord->GetCellVolume(k, j, i);
              }
            }
          }
        }
      }
    }
    // Get total volume for each SN
    MPI_Allreduce(SNiaProp_l.data(), SNiaProp_g.data(), 6 * Globals::nranks * SN_NUM_MAX_PROC, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    SNiaInfo.clear();
    for (int n = 0; n < Globals::nranks * SN_NUM_MAX_PROC; ++n)
    {
      if (SNiaProp_g(n, 3) > TINY_NUMBER)
      {
        Real sn_radius = SNiaProp_g(n, 3);
        Real sn_sph_vol = 4 * PI * CUBE(sn_radius) / 3;
        SNiaProp_g(n, 5) = std::max(SNiaProp_g(n, 5), sn_sph_vol);
        SNiaProp_g(n, 4) = SNiaProp_g(n, 4) / SNiaProp_g(n, 5);
        SNiaInfo.add(SNiaProp_g(n, 0), SNiaProp_g(n, 1), SNiaProp_g(n, 2), SNiaProp_g(n, 3), SNiaProp_g(n, 4), SNiaProp_g(n, 5));
      }
    }
#endif // SNiaFeedbackFlag

#ifdef SNiiFeedbackFlag
    MPI_Allreduce(SNiiProp_l.data(), SNiiProp_g.data(), 6 * Globals::nranks * SN_NUM_MAX_PROC, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    // Compute SNii volume
    for (int n = 0; n < Globals::nranks * SN_NUM_MAX_PROC; ++n)
    {
      Real sn_x1, sn_x2, sn_x3, sn_radius;
      sn_x3 = SNiiProp_g(n, 0);
      sn_x2 = SNiiProp_g(n, 1);
      sn_x1 = SNiiProp_g(n, 2);
      sn_radius = SNiiProp_g(n, 3);
      if (sn_radius > TINY_NUMBER)
      {
        for (int k = ks; k <= ke; ++k)
        {
          for (int j = js; j <= je; ++j)
          {
            for (int i = is; i <= ie; ++i)
            {
              Real x1 = pcoord->x1v(i), x2 = pcoord->x2v(j), x3 = pcoord->x3v(k);
              Real distance2 = SQR(x1) + SQR(sn_x1) - 2 * x1 * sn_x1 * (std::cos(x2) * std::cos(sn_x2) + std::cos(x3 - sn_x3) * std::sin(x2) * std::sin(sn_x2));
              if (distance2 < SQR(sn_radius))
              {
                SNiiProp_l(n, 5) += pcoord->GetCellVolume(k, j, i);
              }
            }
          }
        }
      }
    }
    // Get total volume for each SN
    MPI_Allreduce(SNiiProp_l.data(), SNiiProp_g.data(), 6 * Globals::nranks * SN_NUM_MAX_PROC, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    SNiiInfo.clear();
    for (int n = 0; n < Globals::nranks * SN_NUM_MAX_PROC; ++n)
    {
      if (SNiiProp_g(n, 3) > TINY_NUMBER)
      {
        Real sn_radius = SNiiProp_g(n, 3);
        Real sn_sph_vol = 4 * PI * CUBE(sn_radius) / 3;
        SNiiProp_g(n, 5) = std::max(SNiiProp_g(n, 5), sn_sph_vol);
        SNiiProp_g(n, 4) = SNiiProp_g(n, 4) / SNiiProp_g(n, 5);
        SNiiInfo.add(SNiiProp_g(n, 0), SNiiProp_g(n, 1), SNiiProp_g(n, 2), SNiiProp_g(n, 3), SNiiProp_g(n, 4), SNiiProp_g(n, 5));
      }
    }
#endif // SNiiFeedbackFlag
#endif // SNiaFeedbackFlag || SNiiFeedbackFlag
  }
  return;
}
//========================================================================================
//! \fn void Mesh::UserWorkInLoop()
//! \brief Function called in main after each cycle.
//========================================================================================
void Mesh::UserWorkInLoop()
{
  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//! \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop()
{
  Real dt = pmy_mesh->dt;
  Real time = pmy_mesh->time;
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real gammaGas = peos->GetGamma();
  Real gamma1 = gammaGas - 1.0;
  const AthenaArray<Real> &prim = phydro->w;
  const AthenaArray<Real> &psr = pscalars->r;
#if defined StarEvolutionFlag || defined StarFormationFlag
  for (int k = ks; k <= ke; ++k)
  {
    for (int j = js; j <= je; ++j)
    {
      for (int i = is; i <= ie; ++i)
      {
        Real &star = ruser_meshblock_data[STAR](k, j, i);
        Real &nova = ruser_meshblock_data[NOVA](k, j, i);
        Real &novaDenDelay = ruser_meshblock_data[NOVA_DELAY](k, j, i);
        const Real &rhoDotSF = ruser_meshblock_data[RHODOTSF](k, j, i);
        const Real &rhoDotEV = ruser_meshblock_data[RHODOTEV](k, j, i);
        Real &newStar = ruser_meshblock_data[NEWSTAR](k, j, i);
        star -= dt * rhoDotEV;
        nova = nova * (1 - dt / tauII) + dt * novaDenDelay / tauII;
        novaDenDelay *= (1 - dt / tauDelay);
        novaDenDelay += massiveStar * dt * rhoDotSF / tauDelay; // The rhoDot of SNII parent stars
        star += dt * rhoDotSF * (1 - massiveStar);
        newStar += dt * rhoDotSF;
      }
    }
  }
#endif // StarEvolutionFlag || StarFormationFlag
  // wind refer to Yuan18
  Real &m_bh = pmy_mesh->ruser_mesh_data[0](M_BH);
  Real &mdisk = pmy_mesh->ruser_mesh_data[0](MDISK);
  Real &mfall = pmy_mesh->ruser_mesh_data[0](MFALL);
  Real &mdot_bh = pmy_mesh->ruser_mesh_data[0](MDOT_BH);
  Real &mdot_bondi = pmy_mesh->ruser_mesh_data[0](MDOT_BONDI);
  Real &mdot_wind = pmy_mesh->ruser_mesh_data[0](MDOT_WIND);
  Real &v_wind = pmy_mesh->ruser_mesh_data[0](V_WIND);
  Real l_edd = 1.26e38 * m_bh / punit->solar_mass_code * punit->erg_s_code;
  Real mdot_edd = l_edd / (0.1 * SQR(punit->speed_of_light_code));
  Real eps_wind;
  int mode_wind;
  Real mdot_jet = 0;
  Real v_jet = 0;
  Real eps_jet = 0;
  Real r_g = 2 * punit->grav_const_code * m_bh / SQR(punit->speed_of_light_code);
  Real tau_disk = r_g / (7.6e8 * alpha_disk * std::pow(nschr_in, -3.5) / punit->code_velocity_cgs); // Kato et al. 2008 eq. 3.61
  mdot_bondi = GetMdot_Bondi(pmy_mesh, pcoord, phydro->w, is, js, je, ks, ke);
  Real mdot_crit = MdotWind_Cold(0.02 * mdot_edd) + 0.02 * mdot_edd;
  if (mdot_bondi > mdot_crit) // cold mode
  {
    if (mdot_bondi > 1.66 * mdot_edd) // super Eddington
    {
      mdot_bh = 0.5874 * std::pow(mdot_bondi / mdot_edd, 1.0593) * mdot_edd;
      mdot_wind = mdot_bondi - mdot_bh;
      v_wind = 0.333153 * std::pow(mdot_bh / mdot_edd, -0.0848) * punit->speed_of_light_code;
      mdisk = std::max(0.0, mdisk - dt * mdot_bondi);
      mode_wind = SUP;
    }
    else
    {
      // Real mdot_inn = mdisk * expFactor(dt, tau_disk);
      Real mdot_inn = mdot_bondi;
      mdot_bh = GetMdot_BH_Cold(mdot_inn);
      mdot_wind = mdot_inn - mdot_bh;
      mdisk = std::max(0.0, mdisk - dt * mdot_inn);
      Real l_bh_cold = 0.1 * SQR(punit->speed_of_light_code) * mdot_bh;
      v_wind = 2.5e4 * std::pow(l_bh_cold / (1.e45 * punit->erg_s_code), 0.4) * punit->km_s_code;
      v_wind = std::min(v_wind, 0.3 * punit->speed_of_light_code);
      mode_wind = SUB;
    }
    eps_wind = SQR(130 * punit->km_s_code) / (gamma1 * gammaGas);
  }
  else if (mdot_bondi > TINY_NUMBER) // hot mode
  {
    Real r_s = 2 * punit->grav_const_code * m_bh / SQR(punit->speed_of_light_code);
    Real r_tr = 3 * r_s * SQR(0.02 * mdot_edd / mdot_bondi);
    // refer to Yuan et al. 2018 eq.19. To make mdot matched with radius, if transition radius is larger than inner boundary
    r_tr = std::min(r_tr, x1min);
    r_tr = std::max(r_tr, 3 * r_s);
    Real mdot_bondi_r_tr_is_r_in = 0.02 * std::sqrt(3. * r_s / x1min) * mdot_edd;
    // if (mdot_bondi > 0.011 * mdot_edd) // bridge between cold and hot mode by Bocheng
    // {
    //   mdot_wind = 1.94583 * std::pow((mdot_bondi / mdot_edd), 1.32325) * mdot_edd;
    //   mdot_bh = mdot_bondi - mdot_wind;
    // }
    if (mdot_bondi > mdot_bondi_r_tr_is_r_in)
    {
      // mdot_bh = std::pow(10, std::log10(x1min / 3 / r_s) / std::log10(mdot_crit / mdot_bondi_r_tr_is_r_in) * std::log10(mdot_bondi / mdot_crit)) * 0.02 * mdot_edd;
      // mdot_wind = mdot_bondi - mdot_bh;
      Real x2 = std::log10(mdot_crit / mdot_edd);
      Real x1 = std::log10(mdot_bondi_r_tr_is_r_in / mdot_edd);
      Real a = (2.85 - 2 * std::log10(x1min / 3 / r_s) / (x2 - x1)) / SQR(x2 - x1);
      Real b = (-1.15 - 3 * (SQR(x2) - SQR(x1)) * a) / (x2 - x1) / 2;
      Real c = 0.85 - 3 * SQR(x2) * a - 2 * x2 * b;
      Real d = std::log10(0.02) - (a * CUBE(x2) + b * SQR(x2) + c * x2);
      Real log_mdot_bondi_edd = std::log10(mdot_bondi / mdot_edd);
      mdot_bh = std::pow(10, a * CUBE(log_mdot_bondi_edd) + b * SQR(log_mdot_bondi_edd) + c * log_mdot_bondi_edd + d) * mdot_edd;
      mdot_wind = mdot_bondi - mdot_bh;
    }
    else
    {
      mdot_bh = mdot_bondi * std::sqrt(3. * r_s / r_tr);
      mdot_wind = mdot_bondi - mdot_bh;
    }
    mdisk = std::max(0.0, mdisk - dt * mdot_bondi);
    v_wind = 0.2 * std::sqrt(punit->grav_const_code * m_bh / r_tr);
    eps_wind = 0.5 / (gamma1 * gammaGas) * punit->grav_const_code * m_bh / (3. * r_tr) * std::pow(x1min / r_tr, -2.0 * (gamma_wind - 1.0));
    mode_wind = HOT;
    // TODO: not sure about the jet
    mdot_jet = 0.5 * mdot_bh;
    v_jet = 0.3 * punit->speed_of_light_code;
    eps_jet = 0;
  }
  else
  {
    mdot_wind = 0;
    v_wind = 0;
    mdot_bh = 0;
    mode_wind = NONE;
  }
  m_bh += dt * mdot_bh;

  Real tau_ff_bh = 0.5 * PI * std::sqrt(CUBE(x1min) / (2. * punit->grav_const_code * m_bh));
#ifdef StellarFeedbackFlag
  Real tau_ff_gal = 0.886 * x1min / sigma0;
#else
  Real tau_ff_gal = tau_ff_bh;
#endif
  Real tau_ff = std::min(tau_ff_bh, tau_ff_gal);
  Real mdot_disk = mfall * expFactor(dt, tau_ff);
  mdisk += dt * mdot_disk;
  mfall -= dt * mdot_disk;
  // Real dt_bh_disk = mdot_disk > TINY_NUMBER ? 0.1 * std::min(mdisk / mdot_disk, tau_disk) : 0.1 * tau_disk;
  mfall += dt * mdot_bondi;
  // Real dt_bh_fall = mdot_disk > TINY_NUMBER ? 0.1 * std::min(mfall / mdot_bondi, tau_ff) : 0.1 * tau_ff;

#ifdef AGNWindFlag
  Outflow wind(pmy_mesh->ruser_mesh_data[1], pmy_mesh->iuser_mesh_data[1]);
  Real time_wind = time + x1min / v_wind;
  // set wind with delay time
  if (v_wind > TINY_NUMBER)
  {
    wind.push_back(mdot_wind, mdot_wind * v_wind, eps_wind, dt, time_wind, mode_wind);
  }
#endif // AGNWindFlag
#ifdef AGNJetFlag
  Outflow jet(pmy_mesh->ruser_mesh_data[2], pmy_mesh->iuser_mesh_data[2]);
  Real time_jet = time + x1min / v_jet;
  // set jet with delay time
  if (v_jet > TINY_NUMBER)
  {
    jet.push_back(mdot_jet, mdot_jet * v_jet, eps_jet, dt, time_jet, JET);
  }
#endif // AGNJetFlag
  Real eps_rad = GetRadEfficiency(mdot_bh / mdot_edd);
  Real &l_bh = pmy_mesh->ruser_mesh_data[0](L_BH);
  l_bh = eps_rad * mdot_bh * SQR(punit->speed_of_light_code);
#ifdef AGNRadFlag
  AthenaArray<Real> &agn_flux = ruser_meshblock_data[AGN_FLUX];
  AthenaArray<Real> &heat_agn = ruser_meshblock_data[HEAT_RATE_AGN];
  AthenaArray<Real> &heat_uvb = ruser_meshblock_data[HEAT_RATE_UVB];
  int size = block_size.nx2 * block_size.nx3;
  Real *buf = new Real[size];
  int nbi = 0;
  for (int k = ks; k <= ke; ++k)
  {
    for (int j = js; j <= je; ++j)
    {
      for (int i = is - 1; i <= ie; ++i)
      {
        if (agn_flux(k, j, i) == 0)
        {
          heat_agn(k, j, i) = 0;
          continue;
        }
        Real flux0 = agn_flux(k, j, i);
        Real rho0 = prim(IDN, k, j, i);
        Real Z0 = psr(METAL, k, j, i) / Constants::solar_metallicity;
        Real temp0 = GetTemp(prim(IPR, k, j, i), rho0, Z0);
        heat_agn(k, j, i) = pcool->AGNHeatingRate(rho0, flux0, redshift, Z0, temp0, heat_uvb(k, j, i));
      }
    }
  }
  if (pcoord->x1v(is) < x1min + pcoord->dx1f(is))
  {
    for (int k = ks; k <= ke; k++)
    {
      for (int j = js; j <= je; j++)
      {
        agn_flux(k, j, is - 1) = l_bh / (4 * PI * SQR(pcoord->x1v(is - 1)));
      }
    }
  }
#ifdef MPI_PARALLEL
  else
  {
    const NeighborBlock &nb = pbval->neighbor[nbi];
    MPI_Recv(buf, size, MPI_ATHENA_REAL, nb.snb.rank, 1001, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int p = 0;
    BufferUtility::UnpackData(buf, agn_flux, is - 1, is - 1, js, je, ks, ke, p);
    nbi++;
  }
#endif // MPI_PARALLEL
  for (int k = ks; k <= ke; ++k)
  {
    for (int j = js; j <= je; ++j)
    {
      for (int i = is; i <= ie; ++i)
      {
        if (agn_flux(k, j, i - 1) == 0)
        {
          agn_flux(k, j, i) = 0;
          continue;
        }
        agn_flux(k, j, i) = (agn_flux(k, j, i - 1) * SQR(pcoord->x1v(i - 1)) - (heat_agn(k, j, i - 1) + heat_agn(k, j, i)) / 2 * SQR(pcoord->x1f(i)) * pcoord->dx1v(i - 1)) / SQR(pcoord->x1v(i));
        agn_flux(k, j, i) = std::max(agn_flux(k, j, i), 0.0);
      }
    }
  }
#ifdef MPI_PARALLEL
  if (pcoord->x1v(ie) < x1max - pcoord->dx1f(ie)) // most outer boundary need not pass data.
  {
    const NeighborBlock &nb = pbval->neighbor[nbi];
    int p = 0;
    BufferUtility::PackData(agn_flux, buf, ie, ie, js, je, ks, ke, p);
    MPI_Send(buf, size, MPI_ATHENA_REAL, nb.snb.rank, 1001, MPI_COMM_WORLD);
  }
#endif // MPI_PARALLEL
  delete[] buf;
#endif // AGNRadFlag

#ifdef CheckNaN
  std::stringstream msg;
  msg.precision(3);
  msg.setf(std::ios::fixed);
  for (int k = ks; k <= ke; k++)
  {
    for (int j = js; j <= je; j++)
    {
      for (int i = is; i <= ie; i++)
      {
        std::stringstream error_msg;
        error_msg.precision(3);
        error_msg.setf(std::ios::scientific);

        Real x1 = pcoord->x1v(i) / punit->kpc_code;
        Real x2 = pcoord->x2v(j) * 180. / PI;
        Real x3 = pcoord->x3v(k) * 180. / PI;

        Real vel_mag = std::sqrt(SQR(phydro->w(IVX, k, j, i)) + SQR(phydro->w(IVY, k, j, i)) + SQR(phydro->w(IVZ, k, j, i)));
        if (vel_mag > punit->speed_of_light_code)
        {
          error_msg << std::scientific << "  velocity magnitude    = " << vel_mag << " > c, where v1 = " << phydro->w(IVX, k, j, i) << ", v2 = " << phydro->w(IVY, k, j, i) << ", v3 = " << phydro->w(IVZ, k, j, i) << "\n";
        }

        if (CheckNaN >> 0 & 1)
        {
          const PrimIndex CheckPrimIndex[] = {IVX, IVY, IVZ, IPR};
          const char *PrimIndexString[] = {"IDN", "IVX", "IVY", "IVZ", "IPR"};
          for (auto iter : CheckPrimIndex)
          {
            Real value = phydro->w(iter, k, j, i);
            if (!CheckCondition(iter, value))
            {
              error_msg << "  in PrimIndex          = " << PrimIndexString[iter] << "  found value " << value << "\n";
            }
          }
        }

        if (CheckNaN >> 1 & 1)
        {
          const ConsIndex CheckConsIndex[] = {IDN, IM1, IM2, IM3, IEN};
          const char *ConsIndexString[] = {"IDN", "IM1", "IM2", "IM3", "IEN"};
          for (auto iter : CheckConsIndex)
          {
            Real value = phydro->u(iter, k, j, i);
            if (!CheckCondition(iter, value))
            {
              error_msg << "  in ConsIndex          = " << ConsIndexString[iter] << "  found value " << value << "\n";
            }
          }
        }

        if (CheckNaN >> 2 & 1)
        {
          const ScalarIndex CheckScalarIndex[] = {ISM, CGM, SW, SNIA, SNII, AGNWC, AGNWH, AGNJ};
          const char *ScalarIndexString[] = {"ISM", "CGM", "SW", "SNIA", "SNII", "AGNWC", "AGNWH", "AGNJ"};
          if (NSCALARS > 0)
          {
            for (auto iter : CheckScalarIndex)
            {
              Real value = pscalars->s(iter, k, j, i);
              if (!CheckCondition(iter, value))
              {
                error_msg << "  in ScalarIndex        = " << ScalarIndexString[iter] << "  found value " << value << "\n";
              }
            }
          }
        }

        if (CheckNaN >> 3 & 1)
        {
          const UserMeshBlockIndex CheckUserMeshBlockIndex[] = {STAR, NEWSTAR, NOVA, NOVA_DELAY, AGN_FLUX, COOL_RATE, HEAT_RATE_AGN, HEAT_RATE_UVB, RHODOTSF, RHODOTEV};
          const char *UserMeshBlockIndexString[] = {"STAR", "NEWSTAR", "NOVA", "NOVE_DELAY", "AGN_FLUX", "COOL_RATE", "HEAT_RATE_AGN", "HEAT_RATE_UVB", "RHODOTSF", "RHODOTEV"};
          for (auto iter : CheckUserMeshBlockIndex)
          {
            Real value = ruser_meshblock_data[iter](k, j, i);
            if (!CheckCondition(iter, value))
            {
              error_msg << "  in UserMeshBlockIndex = " << UserMeshBlockIndexString[iter] << "  found value " << value << "\n";
            }
          }
        }

        if (CheckNaN >> 4 & 1)
        {
          if (MAGNETIC_FIELDS_ENABLED)
          {
            Real div = 0, b = 0, ratio = 0, dr = 0;
            div += pcoord->GetFace1Area(k, j, i + 1) * pfield->b.x1f(k, j, i + 1) - pcoord->GetFace1Area(k, j, i) * pfield->b.x1f(k, j, i);
            div += pcoord->GetFace2Area(k, j + 1, i) * pfield->b.x2f(k, j + 1, i) - pcoord->GetFace2Area(k, j, i) * pfield->b.x2f(k, j, i);
            div += pcoord->GetFace3Area(k + 1, j, i) * pfield->b.x3f(k + 1, j, i) - pcoord->GetFace3Area(k, j, i) * pfield->b.x3f(k, j, i);
            div /= pcoord->GetCellVolume(k, j, i);
            b = std::sqrt(SQR(pfield->bcc(IB1, k, j, i)) + SQR(pfield->bcc(IB2, k, j, i)) + SQR(pfield->bcc(IB3, k, j, i)));
            dr = std::sqrt(3. / (1 / SQR(pcoord->dx1f(i)) + 1 / SQR(pcoord->x1v(i) * pcoord->dx2f(j)) + 1 / SQR(pcoord->x1v(i) * std::sin(pcoord->x2v(j)) * pcoord->dx3f(k))));
            ratio = std::abs(div) / (b / dr);
            ruser_meshblock_data[DIVB](k, j, i) = ratio;
            if (ratio > 1e-2)
            {
              error_msg << "  found div B = " << div << ", B/dr = " << (b / dr) << ", ratio = " << ratio << "\n";
            }
          }
        }
        if (error_msg.rdbuf()->in_avail() > 0)
        {
          msg << "  in rank id            = " << Globals::my_rank << "\n"
              << "  in cell (k,j,i)       = (" << k << ", " << j << ", " << i << ")\n"
              << "          (r,theta,phi) = (" << x1 << " kpc, " << x2 << " deg, " << x3 << " deg)\n"
              << error_msg.rdbuf() << "\n";
        }
      }
    }
  }
  if (msg.rdbuf()->in_avail() > 0)
  {
    std::stringstream msg_all;
    msg_all << "\n### INVALID VALUE ###\n"
            << msg.rdbuf()
            << "### PROGRAM ABORT ###\n";
    ATHENA_ERROR(msg_all);
  }
#endif // CheckNaN
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Function called after main loop is finished for user-defined work.
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  // do nothing
  return;
}

//========================================================================================
//! \fn void yuan18Source(MeshBlock *pmb, const Real time, const Real dt,
//!                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//!                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//!                  AthenaArray<Real> &cons_scalar)
//! \brief Add mass and energy return from star formation, star evolution, SN Ia, SN II.
//========================================================================================
void yuan18Source(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar)
{
  Real te_time = std::exp(0.17266 * time / punit->giga_yr_code);
  redshift = 1.22156 / std::pow(0.207671 * te_time + 0.199781 / te_time - 0.407376, 0.333333) - 1; // H0 = 67.74 km/s/Mpc, Omega_m = 0.3089, Omega_Lambda = 0.6911

  Real gammaGas = pmb->peos->GetGamma();
  Real gamma1 = gammaGas - 1.0;
  Coordinates *pcoord = pmb->pcoord;
  AthenaArray<Real> &psr = pmb->pscalars->r;

#ifdef StarEvolutionFlag
  // Stellar mass loss
  if (time <= 0)
  {
    // Abort if the current time is not greater than 0
    std::stringstream msg;
    msg << "Error: Current time should be greater than 0.\n";
    ATHENA_ERROR(msg);
  }
  Real logTime_yr = std::log10(time / punit->yr_code);
  // Eq 4.9, but it is not correct. time should be in years. Refer to Renzini and Buzzoni 1986.
  Real logmTO_sun = 0.0558 * SQR(logTime_yr) - 1.338 * logTime_yr + 7.764;
  Real mTO_sun = std::pow(10, logmTO_sun);
  Real mDotTO = std::abs((mTO_sun * punit->solar_mass_code / time) * (-1.338 + 2 * 0.0558 * logTime_yr));

  Real mRem_sun, mGas_sun;
  if (mTO_sun >= 9)
  {
    mRem_sun = 1.4; // Neutron star
  }
  else
  {
    mRem_sun = 0.503 + 0.055 * mTO_sun;
  }
  mGas_sun = mTO_sun - mRem_sun;

#endif // StarEvolutionFlag

  for (int k = pmb->ks; k <= pmb->ke; ++k)
  {
    for (int j = pmb->js; j <= pmb->je; ++j)
    {
      for (int i = pmb->is; i <= pmb->ie; ++i)
      {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real vol = pcoord->GetCellVolume(k, j, i);
        Real zmetal = psr(METAL, k, j, i) / Constants::solar_metallicity;
        // Take density and temperature at time step n from cons, not from prim
        // Because we do not need intermediate step to calculate the cooling function
        Real &rho = cons(IDN, k, j, i);
        Real &m1 = cons(IM1, k, j, i);
        Real &m2 = cons(IM2, k, j, i);
        Real &m3 = cons(IM3, k, j, i);

        if (rho < pmb->peos->GetDensityFloor())
        {
          printf("WARNING: rho = %e < 0 at (k, j, i) = (%d, %d, %d) in rank %d = (%e kpc, %e deg, %e deg)\n", rho, k, j, i, Globals::my_rank, x1 * punit->kpc_code, x2 * 180. / PI, x3 * 180. / PI);
          Real ek = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
          Real v1 = m1 / rho;
          Real v2 = m2 / rho;
          Real v3 = m3 / rho;
          rho = pmb->peos->GetDensityFloor();
          m1 = v1 * rho;
          m2 = v2 * rho;
          m3 = v3 * rho;
          for (int s = 0; s < NSCALARS; ++s)
          {
            cons_scalar(s, k, j, i) = prim_scalar(s, k, j, i) * rho;
          }
          Real ek_correct = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
          cons(IEN, k, j, i) += ek_correct - ek;
        }
        // Gravity
        Real acc_tot = GetTotalAcceleration(x1);
        m1 += dt * rho * acc_tot;
        if (NON_BAROTROPIC_EOS)
        {
          cons(IEN, k, j, i) += dt * rho * acc_tot * prim(IVX, k, j, i);
        }

        // DO NOT modify name of "ekin_new, eint_new, eint, ekin" to keep the following modules working separately as well
        Real ekin_new = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
        Real eint_new = cons(IEN, k, j, i) - ekin_new;
        Real ekin, eint;
        if (MAGNETIC_FIELDS_ENABLED)
        {
          eint_new -= 0.5 * (SQR(bcc(IB1, k, j, i)) + SQR(bcc(IB2, k, j, i)) + SQR(bcc(IB3, k, j, i)));
        }

#ifdef StellarFeedbackFlag
        Real rhoSNia(0), rhoSNii(0), ergSNia(0), ergSNii(0), ergEV(0), rhoEV(0);
        eint = eint_new;
        ekin = ekin_new;

#ifdef SNiaFeedbackFlag
        for (int n = 0; n < SNiaInfo.num; ++n)
        {
          Real sn_x1, sn_x2, sn_x3, sn_radius;
          sn_x3 = SNiaInfo.info[n].x3;
          sn_x2 = SNiaInfo.info[n].x2;
          sn_x1 = SNiaInfo.info[n].x1;
          sn_radius = SNiaInfo.info[n].radius;
          Real distance2 = SQR(x1) + SQR(sn_x1) - 2 * x1 * sn_x1 * (std::cos(x2) * std::cos(sn_x2) + std::cos(x3 - sn_x3) * std::sin(x2) * std::sin(sn_x2));
          if (distance2 < SQR(sn_radius))
          {
            rhoSNia += SNiaInfo.info[n].rho;
          }
        }
        rhoSNia *= dt / pmb->pmy_mesh->dt;
        ergSNia = rhoSNia * eps_snia;
#ifdef debug
        if (rhoSNia > TINY_NUMBER)
          printf("SNia at (k, j, i) = (%d, %d, %d) in rank %d, (r, theta, phi ) = (%e kpc, %e deg, %e deg), rho ratio = %e, erg ratio = %e\n", k, j, i, Globals::my_rank, x1, x2 * 180. / PI, x3 * 180. / PI, rhoSNia / rho, ergSNia / cons(IEN, k, j, i));
#endif // debug
        if (NSCALARS > 0)
        {
          cons_scalar(SNIA, k, j, i) += rhoSNia;
#ifdef MetalYieldsFlag
          cons_scalar(METAL, k, j, i) += rhoSNia; // SNia composed almost completely of heavy elements
#else
          cons_scalar(METAL, k, j, i) += prim_scalar(METAL, k, j, i) * rhoSNia;
#endif // MetalYieldsFlag
        }
#endif // SNiaFeedbackFlag

#ifdef SNiiFeedbackFlag
        for (int n = 0; n < SNiiInfo.num; ++n)
        {
          Real sn_x1, sn_x2, sn_x3, sn_radius;
          sn_x3 = SNiiInfo.info[n].x3;
          sn_x2 = SNiiInfo.info[n].x2;
          sn_x1 = SNiiInfo.info[n].x1;
          sn_radius = SNiiInfo.info[n].radius;
          Real distance2 = SQR(x1) + SQR(sn_x1) - 2 * x1 * sn_x1 * (std::cos(x2) * std::cos(sn_x2) + std::cos(x3 - sn_x3) * std::sin(x2) * std::sin(sn_x2));
          if (distance2 < SQR(sn_radius))
          {
            rhoSNii += SNiiInfo.info[n].rho;
          }
        }
        rhoSNii *= dt / pmb->pmy_mesh->dt;
        ergSNii = rhoSNii * eps_snii;
#ifdef debug
        if (rhoSNii > TINY_NUMBER)
          printf("SNii at (k, j, i) = (%d, %d, %d) in rank %d, (r, theta, phi ) = (%e kpc, %e deg, %e deg), rho ratio = %e, erg ratio = %e\n", k, j, i, Globals::my_rank, x1, x2 * 180. / PI, x3 * 180. / PI, rhoSNii / rho, ergSNii / cons(IEN, k, j, i));
#endif // debug
        if (NSCALARS > 0)
        {
          cons_scalar(SNII, k, j, i) += rhoSNii;
#ifdef MetalYieldsFlag
          Real ZrhoSNii(0);
          // TODO: consider the metallicities of the SNii progenitors separately--z_snii
          ZrhoSNii = 1.02 * (1.9134 + 0.0479 * std::max(GetStarMetal(x1), 1.65)) / 10.5 * rhoSNii; // FIRE-2 simulation Hopkins et al. (2018b)
          cons_scalar(METAL, k, j, i) += ZrhoSNii;
#else
          cons_scalar(METAL, k, j, i) += prim_scalar(METAL, k, j, i) * rhoSNii;
#endif // MetalYieldsFlag
        }
#endif // SNiiFeedbackFlag

#if defined StarEvolutionFlag || defined StarFormationFlag
        const Real &star = pmb->ruser_meshblock_data[STAR](k, j, i);
#endif // StarEvolutionFlag || StarFormationFlag
#ifdef StarEvolutionFlag
        Real starMass = star * vol;
        Real Lb = starMass / (5.81 * punit->solar_mass_code); // In Lsun
        Real imfTO = 2.9 * Lb * std::pow(mTO_sun, -2.35);     // Salpeter law & C;D;P;R ApJ 376:380, 1991

        Real rhoDotEV = imfTO * mDotTO * mGas_sun / vol;
        pmb->ruser_meshblock_data[RHODOTEV](k, j, i) = rhoDotEV;
        // End Mass loss
#ifndef GravityFlag
        Real rRe = 0.7447 * x1 / rEff;                                                                                                  // Ciotti et al. MNRAS 393:491, 2009, Eq(13)
        Real sigmaR2 = SQR(sigma0) * (6 * SQR(rRe * (1 + rRe)) * std::log((1 + rRe) / rRe) + (1 + rRe) * (1 - 3 * rRe - 6 * SQR(rRe))); // + G * m_bh / (3 * x1) if XPTMASS
                                                                                                                                        // Yoon et al. ApJ 864:6, 2018, Eq(7)
#else
        Real sigmaR2 = 0;
        // TODO follow Ciotti et al. ApJ 933:154, 2022, Eq(A2-A10)
#endif // GravityFlag
        rhoEV = dt * rhoDotEV;
        ergEV = rhoEV * sigmaR2 / gamma1;

        if (NSCALARS > 0)
        {
          cons_scalar(SW, k, j, i) += rhoEV;
#ifdef MetalYieldsFlag
          // TODO: consider the metallicities of the stars as a function of radius. and its evolution.
          cons_scalar(METAL, k, j, i) += rhoEV * GetStarEvoMetal(mTO_sun, GetStarMetal(x1));
#else
          cons_scalar(METAL, k, j, i) += prim_scalar(METAL, k, j, i) * rhoEV;
#endif // MetalYieldsFlag
        }
#endif // StarEvolutionFlag

        eint_new += ergEV + ergSNia + ergSNii;
        rho += rhoEV + rhoSNia + rhoSNii;
        ekin_new = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
        cons(IEN, k, j, i) += ekin_new - ekin + (eint_new - eint);
#endif // StellarFeedbackFlag

        eint = eint_new;
        ekin = ekin_new;
        Real p = eint * gamma1;
        Real temp = GetTemp(p, rho, zmetal);
#if (COOLING_ENABLED == 2)
        Real &cool = pmb->ruser_meshblock_data[COOL_RATE](k, j, i);
        /*if (temp < 0)
        {
          printf("WARNING: temp = %e < 0, press = %e, rho = %e, at (k, j, i) = (%d, %d, %d) in rank %d = (%e kpc, %e deg, %e deg)\n", temp, p, rho, k, j, i, Globals::my_rank, x1 * punit->kpc_code, x2 * 180. / PI, x3 * 180. / PI);
        }*/
#ifdef AGNRadFlag
        Real &agn_flux = pmb->ruser_meshblock_data[AGN_FLUX](k, j, i);
        Real &heat_agn = pmb->ruser_meshblock_data[HEAT_RATE_AGN](k, j, i);

        temp = pmb->pcool->CoolingSrc(rho, agn_flux, redshift, zmetal, temp, dt, cool);
        // Real tcool = eint / cool;
        // user_timestep = std::min(user_timestep, eps_tcool * tcool);
        // rad force in x1 direction: H/c for absorption, flux*rho*kappa_es/c for scattering
        m1 += (heat_agn + agn_flux * rho * kappa_es) / punit->speed_of_light_code * dt;
#else
        temp = pmb->pcool->CoolingSrc(rho, 0, redshift, zmetal, temp, dt, cool);
#endif // AGNRadFlag
        Real p_new = GetPress(temp, rho, zmetal);
        eint_new = p_new / gamma1;
        ekin_new = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
        cons(IEN, k, j, i) += (eint_new - eint) + (ekin_new - ekin);
#endif // COOLING_ENABLED

#ifdef StarFormationFlag
        Real rhoDotSF(0);
        if (temp < temp_sf && rho / GetMeanMolWeight(zmetal) > nrho_sf)
        {
          Real massEnclosed = -GetTotalAcceleration(x1) * SQR(x1) / punit->grav_const_code;
          Real tauRot = 2 * PI * x1 / std::sqrt(punit->grav_const_code * massEnclosed / x1);
          Real tauCollGas = std::sqrt(3 * PI / (32 * punit->grav_const_code * rho));
#if COOLING_ENABLED
          Real tauCool = eint / cool; // eint is corresponding to cool rate.
          if (temp <= pmb->pcool->tfloor)
          {
            tauCool = std::numeric_limits<Real>::min();
          }
#else
          Real tauCool = std::numeric_limits<Real>::min();
#endif // COOLING_ENABLED
       // if temp_new <= tfloor, the cooling does not work, so the tauCool have no meaning, and should be set to a very small value.
          Real tauStar = std::max(tauCool, std::min(tauCollGas, tauRot));
          rhoDotSF = etaSF * rho * expFactor(dt, tauStar);

          if (NSCALARS > 0)
          {
            for (int n = 0; n < NSCALARS; ++n)
            {
              cons_scalar(n, k, j, i) -= dt * rhoDotSF * prim_scalar(n, k, j, i); // Metal will also be subtracted here
            }
            // #ifdef MetalYieldsFlag
            // Real ZDotNSW(0);
            // stellar wind mass loss of new formed star
            // Currently No stellar wind
            // ZDotNSW = (0.0278 + 0.0041 * std::max(zmetal, 1.65)) * 4.763 * (0.01 + zmetal * rhoDotSF / punit->giga_yr_code; //FIRE-2 simulation Hopkins et al. (2018b)
            // NOTE I'm not sure is this ZDotNSW necessary
            // cons_scalar(METAL, k, j, i) += dt * ZDotNSW;
            // #endif // MetalYieldsFlag
          }
          eint = eint_new;
          ekin = ekin_new;
          eint_new -= dt * rhoDotSF * eint / rho;
          m1 -= dt * rhoDotSF * m1 / rho;
          m2 -= dt * rhoDotSF * m2 / rho;
          m3 -= dt * rhoDotSF * m3 / rho;
          rho -= dt * rhoDotSF;
          ekin_new = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
          cons(IEN, k, j, i) += ekin_new - ekin + (eint_new - eint);
        }
        pmb->ruser_meshblock_data[RHODOTSF](k, j, i) = rhoDotSF;
#endif // StarFormationFlag
        bool superluminal = false;
        if (std::abs(m1 / rho) > vmax)
        {
          m1 = (m1 > 0 ? 1 : -1) * vmax * rho;
          superluminal = true;
        }
        if (std::abs(m2 / rho) > vmax)
        {
          m2 = (m2 > 0 ? 1 : -1) * vmax * rho;
          superluminal = true;
        }
        if (std::abs(m3 / rho) > vmax)
        {
          m3 = (m3 > 0 ? 1 : -1) * vmax * rho;
          superluminal = true;
        }
        if (superluminal)
        {
          ekin = ekin_new;
          ekin_new = 0.5 * (SQR(m1) + SQR(m2) + SQR(m3)) / rho;
          cons(IEN, k, j, i) += ekin_new - ekin;
        }
      }
    }
  }
  return;
}

#ifdef CheckNaN
bool CheckCondition(PrimIndex attribute, Real value)
{
  switch (attribute)
  {
  case IPR:
    return value > TINY_NUMBER;
  default:
    return !std::isnan(value);
  }
}

bool CheckCondition(ConsIndex attribute, Real value)
{
  switch (attribute)
  {
  case IDN:
  case IEN:
    return value > TINY_NUMBER;
  default:
    return !std::isnan(value);
  }
}

bool CheckCondition(ScalarIndex attribute, Real value)
{
  switch (attribute)
  {
  default:
    return !std::isnan(value) && value >= 0 && value <= 1;
  }
}
bool CheckCondition(UserMeshBlockIndex attribute, Real value)
{
  switch (attribute)
  {
  case AGN_FLUX:
    return value >= 0;
  default:
    return !std::isnan(value);
  }
}
#endif

void InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
             Real time, Real dt,
             int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  // Ref to TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int stage)
  // only primary variables are used to set boundary condition.
  AthenaArray<Real> &psr = pmb->pscalars->r;
  Real gamma1 = pmb->peos->GetGamma() - 1.0;
#ifdef AGNWindFlag
  Real mdot_wind_now = pmb->pmy_mesh->ruser_mesh_data[0](MDOT_WIND_NOW);
  Real pdot_wind_now = pmb->pmy_mesh->ruser_mesh_data[0](PDOT_WIND_NOW);
  Real eps_wind_now = pmb->pmy_mesh->ruser_mesh_data[0](EPS_WIND_NOW);
  int mode_wind_now = pmb->pmy_mesh->iuser_mesh_data[0](MODE_WIND_NOW);
  Real v_wind_now = pdot_wind_now / mdot_wind_now;
#endif // AGNWindFlag
#ifdef AGNJetFlag
  Real mdot_jet_now = pmb->pmy_mesh->ruser_mesh_data[0](MDOT_JET_NOW);
  Real pdot_jet_now = pmb->pmy_mesh->ruser_mesh_data[0](PDOT_JET_NOW);
  Real eps_jet_now = pmb->pmy_mesh->ruser_mesh_data[0](EPS_JET_NOW);
  int mode_jet_now = pmb->pmy_mesh->iuser_mesh_data[0](MODE_JET_NOW);
  Real v_jet_now = pdot_jet_now / mdot_jet_now;
#endif // AGNJetFlag
  for (int k = ks; k <= ke; ++k)
  {
    for (int j = js; j <= je; ++j)
    {
#if defined(AGNWindFlag) || defined(AGNJetFlag)
      Real area = pco->GetFace1Area(k, j, is);
      int mode = NONE;
      Real frac = 0;
      Real den_outflow = 0, pr_outflow = 0, v_outflow = 0;
      Real jet_ratio = 0;
#endif // AGNWindFlag || AGNJetFlag
#ifdef AGNWindFlag
      if (mode_wind_now == HOT)
      {
        if ((pco->x2v(j) >= ang1_hot && pco->x2v(j) <= ang2_hot) || (pco->x2v(j) <= (PI - ang1_hot) && pco->x2v(j) >= (PI - ang2_hot)))
        {
          frac = area / total_area_hot;
          mode = HOT;
        }
      }
      else if (mode_wind_now == SUB)
      {
        if (ang2_sub > TINY_NUMBER)
        {
          if ((pco->x2v(j) >= ang1_sub && pco->x2v(j) <= ang2_sub) || (pco->x2v(j) <= (PI - ang1_sub) && pco->x2v(j) >= (PI - ang2_sub)))
          {
            frac = area / total_area_sub;
            mode = SUB;
          }
        }
        else
        {
          frac = (CUBE(std::cos(pco->x2f(j))) - CUBE(std::cos(pco->x2f(j + 1)))) * pco->dx3f(k) / PI / 4;
          mode = SUB;
        }
      }
      else if (mode_wind_now == SUP)
      {
        if ((pco->x2v(j) >= ang1_sup && pco->x2v(j) <= ang2_sup) || (pco->x2v(j) <= (PI - ang1_sup) && pco->x2v(j) >= (PI - ang2_sup)))
        {
          frac = area / total_area_sup;
          mode = SUP;
        }
      }
      if (mode)
      {
        den_outflow = mdot_wind_now / v_wind_now / area * frac;
        pr_outflow = gamma1 * eps_wind_now * den_outflow;
        v_outflow = v_wind_now;
      }
#endif // AGNWindFlag
#ifdef AGNJetFlag
      if (mode_jet_now == JET)
      {
        if ((pco->x2v(j) >= ang1_jet && pco->x2v(j) <= ang2_jet) || (pco->x2v(j) <= (PI - ang1_jet) && pco->x2v(j) >= (PI - ang2_jet)))
        {
          if (mode == NONE)
          {
            frac = area / total_area_jet;
            den_outflow = mdot_jet_now / v_jet_now / area * frac;
            pr_outflow = gamma1 * eps_jet_now * den_outflow;
            v_outflow = v_jet_now;
            mode = JET;
          }
          else
          {
            frac = area / total_area_jet;
            Real p_wind = den_outflow * v_outflow;
            Real den_jet = mdot_jet_now / v_jet_now / area * frac;
            Real p_jet = den_jet * v_jet_now;
            den_outflow += den_jet;
            jet_ratio = den_jet / den_outflow;
            pr_outflow = gamma1 * eps_jet_now * den_outflow;
            v_outflow = (p_wind + p_jet) / den_outflow;
            mode += JET;
          }
        }
      }
#endif // AGNJetFlag
#if defined(AGNWindFlag) || defined(AGNJetFlag)
      if (mode)
      {
        for (int i = 1; i <= ngh; i++)
        {
          prim(IVX, k, j, is - i) = v_outflow;
          prim(IDN, k, j, is - i) = den_outflow;
          prim(IPR, k, j, is - i) = pr_outflow > TINY_NUMBER ? pr_outflow : prim(IPR, k, j, is);
          prim(IVY, k, j, is - i) = 0;
          prim(IVZ, k, j, is - i) = 0;
          if (NSCALARS > 0)
          {
            for (int s = 0; s < NSCALARS; s++)
            {
              psr(s, k, j, is - i) = 0.0;
            }
            psr(METAL, k, j, is - i) = z_wind * Constants::solar_metallicity;

            if (mode == HOT)
            {
              psr(AGNWH, k, j, is - i) = 1.0;
            }
            else if (mode == SUB || mode == SUP)
            {
              psr(AGNWC, k, j, is - i) = 1.0;
            }
            else if (mode == JET)
            {
              psr(AGNJ, k, j, is - i) = 1.0;
            }
            else if (mode == JET + SUB || mode == JET + SUP) // Jet conflicts with wind only when jet in hot mode catches up with wind in sub/supersonic mode
            {
              psr(AGNJ, k, j, is - i) = jet_ratio;
              psr(AGNWC, k, j, is - i) = 1.0 - jet_ratio;
            }
          }
        }
      }
      else // no AGN outflow
#endif     // AGNWindFlag || AGNJetFlag
        for (int i = 1; i <= ngh; i++)
        {
          prim(IDN, k, j, is - i) = prim(IDN, k, j, is);
          prim(IVX, k, j, is - i) = std::min(prim(IVX, k, j, is), 0.);
          prim(IVY, k, j, is - i) = prim(IVY, k, j, is);
          prim(IVZ, k, j, is - i) = prim(IVZ, k, j, is);
#ifdef HydrostaticBCFlag
          // integrate pressure to ghost cells to maintain hydrostatic equilibrium
          if (pmb->pcoord->x1f(is - i) <= 0.)
          {
            std::stringstream msg;
            msg << "### FATAL ERROR in function [InnerX1]: "
                << "Non-positive x1f = " << pmb->pcoord->x1f(is - i) << " found at i = " << is - i << "!\n"
                << "Please increase the resolution along r to avoid this!" << std::endl;
            ATHENA_ERROR(msg);
          }
          // This is only first order; higher order might be available for uniform grid and constant gravity
          Real acc_tot = GetTotalAcceleration(pmb->pcoord->x1f(is - i + 1));
          prim(IPR, k, j, is - i) = prim(IPR, k, j, is - i + 1) -
                                    0.5 * (prim(IDN, k, j, is - i) + prim(IDN, k, j, is - i + 1)) *
                                        acc_tot * pco->dx1v(is - i);
#else  // HydrostaticBCFlag not defined
        prim(IPR, k, j, is - i) = prim(IPR, k, j, is);
#endif // HydrostaticBCFlag
          if (NSCALARS > 0)
          {
            for (int s = 0; s < NSCALARS; s++)
            {
              psr(s, k, j, is - i) = psr(s, k, j, is);
            }
          }
        }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED)
  {
    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        for (int i = 1; i <= ngh; ++i)
        {
          b.x1f(k, j, is - i) = b.x1f(k, j, is);
        }
      }
    }
    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je + 1; ++j)
      {
        for (int i = 1; i <= ngh; ++i)
        {
          b.x2f(k, j, is - i) = b.x2f(k, j, is);
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        for (int i = 1; i <= ngh; ++i)
        {
          b.x3f(k, j, is - i) = b.x3f(k, j, is);
        }
      }
    }
  }
}

void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
             Real time, Real dt,
             int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  AthenaArray<Real> &psr = pmb->pscalars->r;
#ifdef CGMInflowFlag
  // Ciotti ApJ 2022 eq 16,19
  Real gamma = pmb->peos->GetGamma();
  Real x1max = pco->x1f(ie + 1);
  Real m_acc_cgm = eps_cgm_mass * Mstar;
  Real mdot_cgm = 2 * m_acc_cgm * time * std::exp(-SQR(time / t0_cgm)) / (SQR(t0_cgm) * (1 - std::exp(-SQR(timespan_cgm / t0_cgm))));
  Real pot_galaxy = PotentialFromGalaxy(x1max);
#ifdef HaloGravityFlag
  Real pot_halo = PotentialFromHalo(x1max);
  Real v_cgm = std::sqrt(-(pot_galaxy + pot_halo) / 2.0);
#else
  Real v_cgm = std::sqrt(-pot_galaxy / 2.0);
#endif // HaloGravityFlag
#endif // CGMInflowFlag

  for (int k = ks; k <= ke; ++k)
  {
    for (int j = js; j <= je; ++j)
    {
#ifdef CGMInflowFlag
      Real frac = 0;
      Real area = pco->GetFace1Area(k, j, is);
      frac = pco->dx3f(k) * ((CUBE(std::cos(pco->x2f(j + 1))) - 3 * std::cos(pco->x2f(j + 1))) - (CUBE(std::cos(pco->x2f(j))) - 3 * std::cos(pco->x2f(j)))) / 8.0 / PI;
      Real den_cgm = mdot_cgm / v_cgm / area * frac;
      Real p_cgm = den_cgm * 1.25 * SQR(v_cgm) / gamma; // 1.25 to avoid supersonic inflow.

#endif // CGMInflowFlag
      for (int i = 1; i <= ngh; i++)
      {
#ifdef CGMInflowFlag
        prim(IDN, k, j, ie + i) = den_cgm;
        prim(IVX, k, j, ie + i) = -v_cgm; // negative for inflow.
        prim(IVY, k, j, ie + i) = prim(IVY, k, j, ie);
        prim(IVZ, k, j, ie + i) = prim(IVZ, k, j, ie);
        prim(IPR, k, j, ie + i) = p_cgm;
        if (NSCALARS > 0)
        {
          for (int s = 0; s < NSCALARS; s++)
          {
            psr(s, k, j, ie + i) = 0.0;
          }
          psr(CGM, k, j, ie + i) = 1.0;
          psr(METAL, k, j, ie + i) = z_cgm * Constants::solar_metallicity;
          ;
        }
#else // CGMInflowFlag not defined
        prim(IDN, k, j, ie + i) = prim(IDN, k, j, ie);
        prim(IVX, k, j, ie + i) = std::max(prim(IVX, k, j, ie), 0.);
        prim(IVY, k, j, ie + i) = prim(IVY, k, j, ie);
        prim(IVZ, k, j, ie + i) = prim(IVZ, k, j, ie);
#ifdef HydrostaticBCFlag
        // integrate pressure to ghost cells to maintain hydrostatic equilibrium
        Real acc_tot = GetTotalAcceleration(pmb->pcoord->x1f(ie + i));
        prim(IPR, k, j, ie + i) = prim(IPR, k, j, ie + i - 1) +
                                  0.5 * (prim(IDN, k, j, ie + i) + prim(IDN, k, j, ie + i - 1)) *
                                      acc_tot * pco->dx1v(ie + i - 1);
        if (prim(IPR, k, j, ie + i) < 0)
        {
          prim(IPR, k, j, ie + i) = prim(IPR, k, j, ie); // In case of extreme situations
        }
#else  // HydrostaticBCFlag not defined
        prim(IPR, k, j, ie + i) = prim(IPR, k, j, ie);
#endif // HydrostaticBCFlag
#endif // CGMInflowFlag
        if (NSCALARS > 0)
        {
          for (int s = 0; s < NSCALARS; s++)
          {
            psr(s, k, j, ie + i) = psr(s, k, j, ie);
          }
        }
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED)
  {
    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        for (int i = 1; i <= ngh; ++i)
        {
          b.x1f(k, j, ie + i + 1) = b.x1f(k, j, ie + 1);
        }
      }
    }
    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je + 1; ++j)
      {
        for (int i = 1; i <= ngh; ++i)
        {
          b.x2f(k, j, ie + i) = b.x2f(k, j, ie);
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        for (int i = 1; i <= ngh; ++i)
        {
          b.x3f(k, j, ie + i) = b.x3f(k, j, ie);
        }
      }
    }
  }
}
void InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
             Real time, Real dt,
             int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  for (int k = ks; k <= ke; ++k)
  {
    for (int i = is; i <= ie; ++i)
    {
      for (int j = 1; j <= ngh; ++j)
      {
        prim(IDN, k, js - j, i) = prim(IDN, k, js, i);
        prim(IVX, k, js - j, i) = prim(IVX, k, js, i);
        prim(IVY, k, js - j, i) = std::min(prim(IVY, k, js, i), 0.);
        prim(IVZ, k, js - j, i) = prim(IVZ, k, js, i);
        prim(IPR, k, js - j, i) = prim(IPR, k, js, i);
        if (NSCALARS > 0)
        {
          for (int s = 0; s < NSCALARS; s++)
          {
            pmb->pscalars->r(s, k, js - j, i) = pmb->pscalars->r(s, k, js, i);
          }
        }
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED)
  {
    for (int k = ks; k <= ke; ++k)
    {
      for (int i = is; i <= ie + 1; ++i)
      {
        for (int j = 1; j <= ngh; ++j)
        {
          b.x1f(k, js - j, i) = b.x1f(k, js, i);
        }
      }
    }
    for (int k = ks; k <= ke; ++k)
    {
      for (int i = is; i <= ie; ++i)
      {
        for (int j = 1; j <= ngh; ++j)
        {
          b.x2f(k, js - j, i) = b.x2f(k, js, i);
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k)
    {
      for (int i = is; i <= ie; ++i)
      {
        for (int j = 1; j <= ngh; ++j)
        {
          b.x3f(k, js - j, i) = b.x3f(k, js, i);
        }
      }
    }
  }
}

void OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
             Real time, Real dt,
             int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  for (int k = ks; k <= ke; ++k)
  {
    for (int i = is; i <= ie; ++i)
    {
      for (int j = 1; j <= ngh; ++j)
      {
        prim(IDN, k, je + j, i) = prim(IDN, k, je, i);
        prim(IVX, k, je + j, i) = prim(IVX, k, je, i);
        prim(IVY, k, je + j, i) = std::max(prim(IVY, k, je, i), 0.);
        prim(IVZ, k, je + j, i) = prim(IVZ, k, je, i);
        prim(IPR, k, je + j, i) = prim(IPR, k, je, i);
        if (NSCALARS > 0)
        {
          for (int s = 0; s < NSCALARS; s++)
          {
            pmb->pscalars->r(s, k, je + j, i) = pmb->pscalars->r(s, k, je, i);
          }
        }
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED)
  {
    for (int k = ks; k <= ke; ++k)
    {
      for (int i = is; i <= ie + 1; ++i)
      {
        for (int j = 1; j <= ngh; ++j)
        {
          b.x1f(k, je + j, i) = b.x1f(k, je, i);
        }
      }
    }
    for (int k = ks; k <= ke; ++k)
    {
      for (int i = is; i <= ie; ++i)
      {
        for (int j = 1; j <= ngh; ++j)
        {
          b.x2f(k, je + j + 1, i) = b.x2f(k, je + 1, i);
        }
      }
    }
    for (int k = ks; k <= ke + 1; ++k)
    {
      for (int i = is; i <= ie; ++i)
      {
        for (int j = 1; j <= ngh; ++j)
        {
          b.x3f(k, je + j, i) = b.x3f(k, je, i);
        }
      }
    }
  }
}

Real UserTimeStep(MeshBlock *pmb)
{
  return user_timestep;
}

Real MeshGeneratorX2(Real x, RegionSize rs)
{
  Real lw, rw;
  Real x2rat = -rs.x2rat;
  if (x2rat == 1.0)
  {
    rw = x, lw = 1.0 - x;
  }
  else
  {
    Real ratn = std::pow(x2rat, rs.nx2 / 2);
    Real rnx = std::pow(x2rat, std::abs((x - 0.5) * rs.nx2));
    lw = 0.5 * (1 - SIGN(x - 0.5) * (rnx - 1) / (ratn - 1));
    rw = 1.0 - lw;
  }
  return rs.x2min * lw + rs.x2max * rw;
}

Real AGNHistoryOutput(MeshBlock *pmb, int iout)
{
  const AthenaArray<Real> &data = pmb->pmy_mesh->ruser_mesh_data[0];
  Real m_bh = data(M_BH);
  Real l_edd = 1.26e38 * m_bh / punit->solar_mass_code * punit->erg_s_code;
  Real mdot_edd = l_edd / (0.1 * SQR(punit->speed_of_light_code));
  switch (iout)
  {
  case 0: // m_bh
    return data(M_BH) / punit->solar_mass_code;
  case 1: // mdot_bondi
    return data(MDOT_BONDI) / mdot_edd;
  case 2: // mdot_bh
    return data(MDOT_BH) / mdot_edd;
  case 3: // mdot_wind
    return data(MDOT_WIND) / mdot_edd;
  case 4: // v_wind
    return data(V_WIND) / punit->speed_of_light_code;
  case 5: // l_bh
#ifdef AGNRadFlag
    return data(L_BH) / l_edd;
#else
    return -data(L_BH) / l_edd; // negative value means the AGN radiation is off
#endif    // AGNRadFlag
  case 6: // mode_now
    return pmb->pmy_mesh->iuser_mesh_data[0](MODE_WIND_NOW) + pmb->pmy_mesh->iuser_mesh_data[0](MODE_JET_NOW);
  case 7: // mdisk
    return data(MDISK) / punit->solar_mass_code;
  case 8: // mfall
    return data(MFALL) / punit->solar_mass_code;
  default:
    return NAN;
  }
}

Real StarHistoryOutput(MeshBlock *pmb, int iout)
{
  switch (iout)
  {
  case agnhst: // star
  {
    Real star = 0;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
    {
      for (int j = pmb->js; j <= pmb->je; ++j)
      {
        for (int i = pmb->is; i <= pmb->ie; ++i)
        {
          Real vol = pmb->pcoord->GetCellVolume(k, j, i);
          star += pmb->ruser_meshblock_data[STAR](k, j, i) * vol;
        }
      }
    }
    return star / punit->solar_mass_code;
  }
  case agnhst + 1: // SFR
  {
    Real sfr = 0;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
    {
      for (int j = pmb->js; j <= pmb->je; ++j)
      {
        for (int i = pmb->is; i <= pmb->ie; ++i)
        {
          Real vol = pmb->pcoord->GetCellVolume(k, j, i);
          sfr += pmb->ruser_meshblock_data[RHODOTSF](k, j, i) * vol;
        }
      }
    }
    return sfr / punit->solar_mass_code * punit->yr_code;
  }
  default:
    return NAN;
  }
}
//========================================================================================
//! \fn Real expFactor(Real dt, Real tau)
//! \brief Generate time step depends on dt and given tau.
//! \details when dt->0, factor->1/tauSF, when dt->inf, factor->1/dt
//! \todo still not sure why need this
//========================================================================================
Real expFactor(Real dt, Real tau)
{
  if (dt < 0.01 * tau)
  {
    return (1 - 0.5 * dt / tau) / tau;
  }
  else
  {
    return (1 - std::exp(-dt / tau)) / dt;
  }
}

// SECTION Potential
Real GetTotalAcceleration(Real r)
{
  Real acc_tot;
  acc_tot = AccelerationFromGalaxy(r);
#ifdef GravityFlag
  Real acc_bh = AccelerationFromBH(r);
  acc_tot += acc_bh;
#ifdef HaloGravityFlag
  Real acc_halo = AccelerationFromHalo(r);
  acc_tot += acc_halo;
#endif // HaloGravityFlag
#endif // GravityFlag
  return acc_tot;
}

Real PotentialFromGalaxy(Real r)
{
#ifdef GravityFlag
  Real s = r / rstar;
  Real pot = G * Mstar * eps_galaxy_mass / (rstar * eps_galaxy_radius) * std::log(s / (eps_galaxy_radius + s));
#else
  Real pot = SQR(vc) * std::log(r / rstar);
#endif // GravityFlag
  return pot;
}

Real AccelerationFromGalaxy(Real r)
{
#ifdef GravityFlag
  Real acc = -G * Mstar * eps_galaxy_mass / (r * (rstar * eps_galaxy_radius + r));
#else
  Real acc = -SQR(vc) / r;
#endif // GravityFlag
  return acc;
}

#ifdef GravityFlag
Real PotentialFromBH(Real r)
{
  Real pot = -G * m_bh / r;
  return pot;
}

Real AccelerationFromBH(Real r)
{
  Real acc = -G * m_bh / SQR(r);
  return acc;
}

#ifdef HaloGravityFlag
Real PotentialFromHalo(Real r)
{
  Real s = r / rstar;
  Real pot = SQR(v_halo) * (0.5 * std::log(1 + SQR(s / eps_halo_radius)) - 1 + eps_halo_radius / s * std::atan(s / eps_halo_radius));
  return pot;
}

Real AccelerationFromHalo(Real r)
{
  Real s = r / rstar;
  Real acc = -SQR(v_halo) * (1 - eps_halo_radius / s * std::atan(s / eps_halo_radius)) / r;
  return acc;
}
#endif // HaloGravityFlag
#endif // GravityFlag
// #SECTION Potential

// SECTION AGN
Real GetMdot_Bondi(Mesh *pm, Coordinates *pco, const AthenaArray<Real> &prim, int is, int js, int je, int ks, int ke)
{
  Real l_mdot = 0;
  Real dm;
  Real area;
  Real mdot = 0;
  if (pco->x1v(is) < pm->mesh_size.x1min + pco->dx1f(is))
  {
    // ibg is the index of the bondi radius(almost) of the disk,
    // mainly used to avoid the unphysical oscillation at the boundary of the domain.
    int ibg = is + 1;
    for (int k = ks; k <= ke; ++k)
    {
      for (int j = js; j <= je; ++j)
      {
        area = pco->GetFace1Area(k, j, ibg);
        dm = -0.5 * area * (prim(IDN, k, j, ibg) * prim(IVX, k, j, ibg) + prim(IDN, k, j, ibg - 1) * prim(IVX, k, j, ibg - 1));
        if (dm > TINY_NUMBER)
        {
          l_mdot += dm;
        }
      }
    }
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(&l_mdot, &mdot, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#else
  mdot = l_mdot;
#endif
  return mdot;
}

Real GetMdot_BH_Cold(Real mdot_in)
{
  Real mdot_bh;
  std::stringstream msg;
  auto func = [mdot_in](Real mdot_bh)
  {
    return MdotWind_Cold(mdot_bh) - (mdot_in - mdot_bh);
  };
  if (Newtonian_Solver(0.5 * mdot_in, 1e-7, func, &mdot_bh) && mdot_bh > TINY_NUMBER && mdot_bh < mdot_in)
  {
    return mdot_bh;
  }
  else
  {
    printf("mdot_bh = %e, mdot_in = %e\n", mdot_bh, mdot_in);
    msg << "### FATAL ERROR in function [yuan18::GetMdot_BH]"
        << std::endl
        << "Mdot_BH could not be calculated";
    ATHENA_ERROR(msg);
  }
}

Real MdotWind_Cold(Real mdot_bh)
{
  Real l_bh_cold = 0.1 * mdot_bh * SQR(punit->speed_of_light_code);
  Real mdot_wind = 0.28 * std::pow(l_bh_cold / (1e45 * punit->erg_s_code), 0.85) * punit->solar_mass_code / punit->yr_code;
  return mdot_wind;
}

int Newtonian_Solver(Real x0, Real eps, const std::function<Real(Real)> &f, Real *root)
// x0 is initial value，f() is the function to be solved
{
  Real x1;
  for (int i = 0; i < 1e5; i++)
  {
    x1 = x0 - f(x0) / ((f(x0 + eps * x0) - f(x0 - eps * x0)) / (2 * eps * x0));
    if (std::abs(x1 - x0) / x0 < eps)
    {
      *root = x1;
      return 1;
    }
    x0 = x1;
  }
  return 0;
}

// Calculate metal source term from star evolution based on bilinear interpolation of
// yields tables from Nomoto et al. (2013)
Real GetStarEvoMetal(Real mTO_sun, Real z_star)
{
  Real Z_enrich_SE(0), Z_SE(0);
  Real Z_star = z_star * Constants::solar_metallicity;

  Real inte1(0), inte2(0);
  int ni1, ni2, nj1, nj2;

  if (mTO_sun <= m_star[0])
  {
    nj1 = 0;
    nj2 = 1;
  }
  else if (mTO_sun >= m_star[9])
  {
    nj1 = 8;
    nj2 = 9;
  }
  else
  {
    for (int j = 0; j < 10; ++j)
    {
      if (mTO_sun < m_star[j])
      {
        nj1 = j - 1;
        nj2 = j;
        break;
      }
    }
  }

  if (Z_star <= Z_s[0])
  {
    ni1 = 0;
    ni2 = 1;
  }
  else if (Z_star >= Z_s[5])
  {
    ni1 = 4;
    ni2 = 5;
  }
  else
  {
    for (int i = 0; i < 6; ++i)
    {
      if (Z_star < Z_s[i])
      {
        ni1 = i - 1;
        ni2 = i;
        break;
      }
    }
  }

  inte1 = (yields[ni2][nj1] - yields[ni1][nj1]) / (Z_s[ni2] - Z_s[ni1]) * (Z_star - Z_s[ni1]) + yields[ni1][nj1];
  inte2 = (yields[ni2][nj2] - yields[ni1][nj2]) / (Z_s[ni2] - Z_s[ni1]) * (Z_star - Z_s[ni1]) + yields[ni1][nj2];

  Z_enrich_SE = (inte2 - inte1) / (m_star[nj2] - m_star[nj1]) * (mTO_sun - m_star[nj1]) + inte1;
  if (Z_enrich_SE < 0.)
  {
    Z_enrich_SE = 0.;
  }

  Z_SE = Z_star + Z_enrich_SE;
  return Z_SE;
}

Real GetStarMetal(Real r)
{
  Real z_star = 0;
  if (r < rEff / 8)
  {
    z_star = z_init_in;
  }
  else
  {
    z_star = z_init_eff * std::pow(r / rEff, z_slope);
  }
  return z_star;
}

Real GetRadEfficiency(Real mdot) // in Eddington
{                                // refer to Yuan et al.2018 eq.25
  Real epsilon;
  if (mdot < 9.4e-5)
  {
    epsilon = 0.12 * std::pow(100 * mdot, 0.59);
  }
  else if (mdot >= 9.4e-5 && mdot < 5.0e-3)
  {
    epsilon = 0.026 * std::pow(100 * mdot, 0.27);
  }
  else if (mdot >= 5.0e-3 && mdot < 6.19e-3) // small modification to avoid the discontinuity
  {
    epsilon = 0.5 * std::pow(100 * mdot, 4.53);
  }
  else if (mdot >= 6.19e-3 && mdot < 1)
  {
    epsilon = 0.057;
  }
  else
  {
    epsilon = 0.1197 * std::pow(100 * mdot, -0.17);
  }
  return eff_em_factor / 0.057 * epsilon;
}
// #SECTION AGN

Real GetMeanMolWeight(Real Z_Zsun)
{
  // TODO add the calculation of mean molecular weight based on metallicity and temperature
  // The following assumes full ionization, with X and Z as hydrogen and metal mass fractions
  Real Z_frac = Constants::solar_metallicity * Z_Zsun;
  Real mu = 1.0 / (2. * H_frac + 3. / 4. * (1. - H_frac - Z_frac) + Z_frac / 2.); // mu ~ 0.61633
  return mu * punit->hydrogen_mass_code;
}

Real GetTemp(Real press, Real rho, Real Z_Zsun)
{
  // TODO calculate temp with metallicity.
  Real mean_mol_weight = GetMeanMolWeight(Z_Zsun);
  Real temp = press * mean_mol_weight / (rho * punit->k_boltzmann_code);
  return temp;
}
Real GetPress(Real temp, Real rho, Real Z_Zsun)
{
  // TODO calculate temp with metallicity.
  Real mean_mol_weight = GetMeanMolWeight(Z_Zsun);
  Real press = rho * punit->k_boltzmann_code * temp / mean_mol_weight;
  return press;
}
// ANCHOR After finishing all work, please check all the comments
