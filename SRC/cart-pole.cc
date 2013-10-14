// Copyright (C) 2012 by Antonio El Khoury.

/**
 * \file SRC/cart-pole.cc
 *
 * \brief cart-pole model in MUSCOD.
 */

#include <cmath>
#include <string>
#include <map>
#include <fstream>
#include <sys/time.h>
#include <vector>

#include "def_usrmod.hpp"

std::ofstream posFile;
std::ofstream velFile;
std::ofstream accFile;
std::ofstream waistFile;
std::ofstream hipFile;
std::ofstream zmpabsFile;
std::ofstream zmpforceFile;
std::ofstream zmpFile;
std::ofstream comFile;
std::ofstream forceFile;
std::ofstream cstrFile;
std::ofstream configFile;
std::ofstream controlFile;
std::ofstream diffStateFile;

/** \brief Objective function (Lagrangian type) */
static void lfcn(double *t, double *xd, double *xa, double *u,
		 double *p, double *lval, double *rwh, long *iwh, long *info);

/** \brief Objective function (Mayer type) */
static void mfcn(double *ts, double *sd, double *sa, double *p,
		 double *pr, double *mval, long *dpnd, long *info);

/** \brief Right hand side of the differential equation */
static void ffcn(double *t, double *xd, double *xa, double *u,
		 double *p, double *rhs, double *rwh, long *iwh, long *info);

/** \brief Constraints at interior shooting nodes */
static void rdfcn_i (double *ts, double *sd, double *sa, double *u,
		     double *p, double *pr, double *res, long *dpnd, long *info);

// \brief Plotter function.
static void mplo ( double *t, double *sd, double *sa, double *u,
		   double *p, double *rwh, long *iwh );

// \brief Final results exporter function.
static void mout (long *imos, long *imsn, double *ts, double *te,
		  double *sd, double *sa,  double *u, double *udot,
		  double *ue, double *uedot, double *p, double *pr,
		  double *ccxd, double *mul_ccxd, double *ares,
		  double *mul_ares,double *rd, double *mul_rd,
		  double *rc, double *mul_rc, double *obj,
		  double *rwh, long *iwh);

long TICKS_PER_SECOND = 1e6;
struct timeval tv_start, tv_stop;

long time_usec = 0;
long inner_loop_time;

long nb_calls = 0;

/** \brief Entry point for the muscod application */
extern "C" void def_model(void);
void def_model(void)
{
  // Initialize output files.
  // posFile.open ("trajectory.pos");
  // velFile.open ("trajectory.vel");
  // accFile.open ("trajectory.acc");
  // waistFile.open ("trajectory.waist");
  // hipFile.open ("trajectory.hip");
  // zmpabsFile.open ("trajectory.zmpabs");
  // zmpforceFile.open ("trajectory.zmpforce");
  // zmpFile.open ("trajectory.zmp");
  // comFile.open ("trajectory.com");
  // forceFile.open ("trajectory.force");
  // cstrFile.open ("trajectory.cstr");
  // configFile.open ("trajectory.config");
  // controlFile.open ("trajectory.control");
  // diffStateFile.open ("trajectory.diffstate");

#define  NMOS   1  /* Number of phases (Model Stages) */
#define  NP     0  /* Number of parameters */
#define  NRC    0  /* Number of coupled constraints */
#define  NRCE   0  /* Number of coupled equality constraints */

#define  NXD    4 /* Number of differential states */
#define  NXA    0 /* Number of algebraic states */
#define  NU     1 /* Number of controls */
#define  NPR    0  /* Number of local parameters */

/* Number of constraints at interior shooting nodes */
#define  NRD_i  0
/* Number of equality constraints at interior shooting nodes */
#define  NRDE_i 0

  /* Define problem dimensions */
  def_mdims(NMOS, NP, NRC, NRCE);
  /* Define the first (and only) phase */
  def_mstage(
	     0,
	     NXD, NXA, NU,
	     NULL, lfcn,
	     0, 0, 0, NULL, ffcn, NULL,
	     NULL, NULL
	     );
  /* Define constraints at all shooting nodes */
  // def_mpc(0, "s", NPR, NRD_i, NRDE_i, rdfcn_i, NULL);
  // def_mpc(0, "i", NPR, NRD_i, NRDE_i, rdfcn_i, NULL);
  // def_mpc(0, "e", NPR, NRD_i, NRDE_i, rdfcn_i, NULL);
  // Define function that will output the results after the
  // opimization is done.
  def_mio (NULL, mout, mplo);

  // std::cout << "NXD " << NXD << std::endl;
  // std::cout << "NU " << NU << std::endl;
}

struct Sample
{
  double t;
  std::vector<double> sd;
  std::vector<double> u;
};

typedef std::map<double, Sample> timeToSample_t;
timeToSample_t timeToSample;

// Cost function parameters.
static const double k_force = 1;
static const double k_height = 1;
static const double k_position = 0.1;
static const double k_velocity = 0.001;

static void lfcn(double *t, double *xd, double *xa, double *u,
		 double *p, double *lval, double *rwh, long *iwh, long *info)
{
  lval[0] = k_force * u[0] * u[0]
    + k_position * xd[0] * xd[0]
    + k_height * (1 - cos (xd[1]))
    + k_velocity * xd[3] * xd[3];
}

// Dynamic system parameters.
static const double g = 9.81;
static const double m = 5.1057;
static const double M = 9.4248;
static const double l = 0.3;
static const double mu_th = 0.1; // damping on \dot{\theta}
static const double gear_ratio = 100.; // Actuator reduction ratio

static void ffcn(double *t, double *xd, double *xa, double *u,
		 double *p, double *rhs, double *rwh, long *iwh, long *info)
{
  rhs[0] = xd[2];
  rhs[1] = xd[3];
  rhs[3] = ((gear_ratio * u[0] - m * l * xd[3] * xd[3] * sin (xd[1]))
	    / (M + m) * cos (xd[1])
	    + g * sin (xd[1]) - mu_th / (m * l) * xd[1])
    / l / (1 - m / (M + m) * cos (xd[1]) * cos (xd[1]));
  rhs[2] = (u[0] + m * l * (rhs[3] * cos (xd[1]) - xd[3] * xd[3] * sin (xd[1])))
    / (M + m);
}

static void mplo ( double *t, double *sd, double *sa, double *u,
		   double *p, double *rwh, long *iwh )
{
  // long conv_achieved_flag;
  // convergenceAchieved (&conv_achieved_flag);
  // if (!conv_achieved_flag)
  //   return;

  // If a new SQP iteration has started, erase data and write
  // again. Apparently the convergenceAchieved function does not work
  // properly.
  if (t[0] == 0)
    {
      timeToSample.clear ();
    }

  // Save all raw data and sample it later.
  Sample sample;
  sample.t = t[0];
  sample.sd.resize (NXD);
  sample.u.resize (NU);
  for (unsigned i = 0; i < NXD; ++i)
    sample.sd[i] = sd[i];
  for (unsigned i = 0; i < NU; ++i)
    sample.u[i] = u[i];

  timeToSample[t[0]] = sample;
}

static void mout (long *imos, long *imsn, double *ts, double *te,
		  double *sd, double *sa,  double *u, double *udot,
		  double *ue, double *uedot, double *p, double *pr,
		  double *ccxd, double *mul_ccxd, double *ares,
		  double *mul_ares, double *rd, double *mul_rd,
		  double *rc, double *mul_rc, double *obj,
		  double *rwh, long *iwh)
{
  // Export Seqplay files. Make sure this function is called just once
  // as all data is treated in on call.
  if (imsn[0] != 0)
    return;

  // Run a fist loop to output raw data in the form:
  // full-configuration
  for (timeToSample_t::iterator it = timeToSample.begin ();
       it != timeToSample.end (); ++it)
    {
      double time = it->first;
      std::vector<double> diffState = it->second.sd;
      std::vector<double> control = it->second.u;

      // Control file
      controlFile << time << " ";
      for (unsigned i = 0; i < NU - 1; ++i)
      	controlFile << control[i] << " ";
      controlFile << control[NU - 1] << "\n";

      // Differential state file
      diffStateFile << time << " ";
      for (unsigned i = 0; i < NXD - 1; ++i)
	diffStateFile << diffState[i] << " ";
      diffStateFile << diffState[NXD - 1] << "\n";
    }
  
  // Get data from map at a sampling period of 5ms.
  double samplingPeriod = 5e-3;
  unsigned sampleId = 0;
  for (timeToSample_t::iterator it = timeToSample.begin ();
       it != timeToSample.end (); ++it)
    {
      double sampleTime = sampleId * samplingPeriod;
      timeToSample_t::iterator sampleIt
	= timeToSample.lower_bound (sampleTime);
      if (sampleIt == timeToSample.end ())
	break;

      // Do a linear interpolation of results if sample time is not
      // found.
      std::vector<double> sampleState = sampleIt->second.sd;
      std::vector<double> sampleControl = sampleIt->second.u;

      if (fabs (sampleIt->first - sampleTime) < 1e-6)
      	{
      	}
      else
      	{
	  if (sampleIt->first > sampleTime)
	    --sampleIt;

      	  timeToSample_t::iterator nextSampleIt = sampleIt;
      	  ++nextSampleIt;
      	  std::vector<double> nextSampleSate = nextSampleIt->second.sd;
      	  std::vector<double> nextSampleControl = nextSampleIt->second.u;

      	  double lambda = (sampleTime - sampleIt->first)
      	    / (nextSampleIt->first - sampleIt->first);
      	  for (unsigned i = 0; i < sampleState.size (); ++i)
      	    sampleState[i] = (1 - lambda) * sampleState[i]
      	      + lambda * nextSampleSate[i];
      	  for (unsigned i = 0; i < sampleControl.size (); ++i)
      	    sampleControl[i] = (1 - lambda) * sampleControl[i]
      	      + lambda * nextSampleControl[i];
      	}

      ++sampleId;
    }

  // Close files.
  posFile.close ();
  velFile.close ();
  accFile.close ();
  waistFile.close ();
  hipFile.close ();
  zmpabsFile.close ();
  zmpFile.close ();
  comFile.close ();
  cstrFile.close ();
  configFile.close ();
  controlFile.close ();
  diffStateFile.close ();
}
