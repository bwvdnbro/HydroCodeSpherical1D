/*******************************************************************************
 * This file is part of HydroCodeSpherical1D
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * HydroCodeSpherical1D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HydroCodeSpherical1D is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with HydroCodeSpherical1D. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "RiemannSolver.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#define RMIN 0.1
#define RMAX 0.475

#define NCELL 100
#define GAMMA 1.001
#define DT 0.00001
#define NSTEP 100000
#define SNAPSTEP 1000

#define NSTEP_RELAX 20000

#define BOUNDARIES_OPEN 0
#define BOUNDARIES_REFLECTIVE 1
#define BOUNDARIES_BONDI 2

#define BOUNDARIES BOUNDARIES_BONDI

#define EOS_IDEAL 0
#define EOS_ISOTHERMAL 1
#define EOS_BONDI 2

#define EOS EOS_BONDI

//#define ISOTHERMAL_U 20.26785
#define ISOTHERMAL_U 31.25

#define POTENTIAL_NONE 0
#define POTENTIAL_POINT_MASS 1

#define POTENTIAL POTENTIAL_POINT_MASS

//#define G 4.30097e-03
#define G 1.

#define MASS_POINT_MASS 1.5

//#define NO_GRADIENTS

/// end of customizable parameters

#define BOXSIZE (RMAX - RMIN)
#define HALF_CELLSIZE (0.5*BOXSIZE/NCELL)
#define ISOTHERMAL_C (ISOTHERMAL_U * (GAMMA - 1.))

class Cell{
public:
  double _rho;
  double _u;
  double _P;

  double _grad_rho;
  double _grad_u;
  double _grad_P;

  double _m;
  double _p;
  double _E;

  double _V;
  double _Vreal;
  double _midpoint;

  double _a;

  Cell *_right_ngb;
};

double density(double r){
  return std::exp(G * MASS_POINT_MASS * (1. / r - 1. / RMIN) / ISOTHERMAL_C);
}

void init(Cell &cell){
//  if(cell._midpoint < 0.25){
//    cell._rho = 1.;
//    cell._P = 1.;
//  } else {
//    cell._rho = 0.125;
//    cell._P = 0.1;
//  }
//  cell._rho = density(cell._midpoint);
//  cell._u = 0.;
//  cell._a = 0.;
  const double r_inv = 0.9 / cell._midpoint;
  cell._rho = std::pow(r_inv, 1.2);
  cell._u = -std::pow(r_inv, 0.8);
  cell._P = 0.03121878122 * cell._rho;
  const double r2 = cell._midpoint * cell._midpoint;
  cell._a = -G * MASS_POINT_MASS / r2;
}

void write_snapshot(unsigned int istep, const Cell cells[NCELL+2]){
  std::stringstream filename;
  filename << "snap_";
  filename.fill('0');
  filename.width(3);
  filename << (istep/SNAPSTEP);
  filename << ".txt";
  std::ofstream ofile(filename.str().c_str());
  for(unsigned int i = 1; i < NCELL+1; ++i){
    ofile << cells[i]._midpoint << "\t" << cells[i]._rho << "\t"
          << cells[i]._u << "\t" << cells[i]._P << "\t" << cells[i]._a << "\n";
  }
  ofile.close();
}

int main(int argc, char **argv){
  std::cout << "ISOTHERMAL_C: " << ISOTHERMAL_C << std::endl;

  Cell cells[NCELL+2];

  for(unsigned int i = 0; i < NCELL+2; ++i){
    cells[i]._midpoint = RMIN + (i - 0.5) * BOXSIZE / NCELL;
    cells[i]._V = BOXSIZE / NCELL;
    const double xmin = RMIN + (i - 1.) * BOXSIZE / NCELL;
    const double xmax = RMIN + i * BOXSIZE / NCELL;
    cells[i]._Vreal = 4. * M_PI * (xmax * xmax * xmax - xmin * xmin * xmin) / 3.;
    if(i < NCELL+1){
      cells[i]._right_ngb = &cells[i+1];
    }
  }

  for(unsigned int i = 1; i < NCELL+1; ++i){
    init(cells[i]);

#if EOS == EOS_ISOTHERMAL || EOS == EOS_BONDI
    cells[i]._P = ISOTHERMAL_C * cells[i]._rho;
#endif

    cells[i]._m = cells[i]._rho * cells[i]._V;
    cells[i]._p = cells[i]._m * cells[i]._u;
    cells[i]._E = cells[i]._P * cells[i]._V / (GAMMA - 1.) +
                  0.5 * cells[i]._u * cells[i]._p;
  }

  RiemannSolver solver(GAMMA);

  for(unsigned int istep = 0; istep < NSTEP; ++istep){

    // add spherical source term
    for(unsigned int i = 1; i < NCELL+1; ++i){
      const double r = cells[i]._midpoint;
      const double Ui[3] = {cells[i]._m / cells[i]._V,
                            cells[i]._p / cells[i]._V,
                            cells[i]._E / cells[i]._V};
      double K1[3];
      K1[0] = -DT * Ui[1] / r;
      K1[1] = -DT * Ui[1] * Ui[1] / Ui[0] / r;
      const double p1 = (GAMMA - 1.) * (Ui[2] - 0.5 * Ui[1] * Ui[1] / Ui[0]);
      K1[2] = -DT * Ui[1] * (Ui[2] + p1) / Ui[0] / r;
      double U[3] = {Ui[0] + K1[0],
                     Ui[1] + K1[1],
                     Ui[2] + K1[2]};
      double K2[3];
      K2[0] = -DT * U[1] / r;
      K2[1] = -DT * U[1] * U[1] / U[0] / r;
      const double p2 = (GAMMA - 1.) * (U[2] - 0.5 * U[1] * U[1] / U[0]);
      K2[2] = -DT * U[1] * (U[2] + p2) / U[0] / r;
      U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);
      U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);
      U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);

      cells[i]._m = U[0] * cells[i]._V;
      cells[i]._p = U[1] * cells[i]._V;
      cells[i]._E = U[2] * cells[i]._V;
    }

    // add gravitational acceleration
#if POTENTIAL == POTENTIAL_POINT_MASS
    for(unsigned int i = 1; i < NCELL+1; ++i){
      const double r = cells[i]._midpoint;
      const double a = -G * MASS_POINT_MASS / (r * r);
      cells[i]._a = a;
      const double m = cells[i]._V * cells[i]._rho;
      cells[i]._p += 0.5 * DT * a * m;
    }
#endif

#if EOS == EOS_BONDI
    double Cion = 5.*std::pow(0.9, 2.4)/3. *
                    (std::pow(0.25, 0.6) - std::pow(0.1, 0.6));
#endif
    for(unsigned int i = 1; i < NCELL+1; ++i){
      cells[i]._rho = cells[i]._m / cells[i]._V;
      cells[i]._u = cells[i]._p / cells[i]._m;
#if EOS == EOS_IDEAL
      cells[i]._P = (GAMMA - 1.) * (cells[i]._E / cells[i]._V -
                              0.5 * cells[i]._rho * cells[i]._u * cells[i]._u);
#elif EOS == EOS_ISOTHERMAL
      cells[i]._P = ISOTHERMAL_C * cells[i]._rho;
#elif EOS == EOS_BONDI
      double sound_speed = ISOTHERMAL_C;
      if(istep > NSTEP_RELAX){
        if(Cion > 0.){
          // assume ionized
          sound_speed *= 10.;
          const double rmin = cells[i]._midpoint - HALF_CELLSIZE;
          const double rmax = cells[i]._midpoint + HALF_CELLSIZE;
          const double Vshell = (rmax * rmax * rmax - rmin * rmin * rmin) / 3.;
          Cion -= Vshell * cells[i]._rho * cells[i]._rho;
        }
      }
      cells[i]._P = sound_speed * cells[i]._rho;
#else
#error "No equation of state selected!"
#endif
    }

    if(istep%SNAPSTEP == 0){
      std::cout << "step " << istep << " of " << NSTEP << std::endl;
      write_snapshot(istep, cells);
    }

#if BOUNDARIES == BOUNDARIES_BONDI
    const double r_inv_low = 0.9 / cells[0]._midpoint;
    cells[0]._rho = std::pow(r_inv_low, 1.2);
    cells[0]._u = -std::pow(r_inv_low, 0.8);
    cells[0]._P = 0.03121878122 * cells[0]._rho;

    const double r_inv_high = 0.9 / cells[NCELL+1]._midpoint;
    cells[NCELL+1]._rho = std::pow(r_inv_high, 1.2);
    cells[NCELL+1]._u = -std::pow(r_inv_high, 0.8);
    cells[NCELL+1]._P = 0.03121878122 * cells[NCELL+1]._rho;
#else
    // the lower boundary is always reflective, as this corresponds to the
    // center of our spherical coordinates
    cells[0]._rho = density(cells[0]._midpoint);
    cells[0]._u = 0.;
    cells[0]._P = ISOTHERMAL_C * cells[0]._rho;
#if BOUNDARIES == BOUNDARIES_OPEN
    cells[NCELL+1]._rho = cells[NCELL]._rho;
    cells[NCELL+1]._u = cells[NCELL]._u;
    cells[NCELL+1]._P = ISOTHERMAL_C * cells[NCELL+1]._rho;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
    cells[NCELL+1]._rho = cells[NCELL]._rho;
    cells[NCELL+1]._u = -cells[NCELL]._u;
    cells[NCELL+1]._P = cells[NCELL]._P;
#else
#error "No boundary conditions chosen!"
#endif
#endif

    for(unsigned int i = 1; i < NCELL+1; ++i){
      double rhomin, rhoplu, umin, uplu, Pmin, Pplu, dmin, dplu;
      rhomin = cells[i]._rho - cells[i-1]._rho;
      rhoplu = cells[i]._rho - cells[i+1]._rho;
      umin = cells[i]._u - cells[i-1]._u;
      uplu = cells[i]._u - cells[i+1]._u;
      Pmin = cells[i]._P - cells[i-1]._P;
      Pplu = cells[i]._P - cells[i+1]._P;
      dmin = cells[i]._midpoint - cells[i-1]._midpoint;
      dplu = cells[i]._midpoint - cells[i+1]._midpoint;
      if(std::abs(rhomin*dplu) < std::abs(rhoplu*dmin)){
        cells[i]._grad_rho = rhomin / dmin;
      } else {
        cells[i]._grad_rho = rhoplu / dplu;
      }
      if(std::abs(umin*dplu) < std::abs(uplu*dmin)){
        cells[i]._grad_u = umin / dmin;
      } else {
        cells[i]._grad_u = uplu / dplu;
      }
      if(std::abs(Pmin*dplu) < std::abs(Pplu*dmin)){
        cells[i]._grad_P = Pmin / dmin;
      } else {
        cells[i]._grad_P = Pplu / dplu;
      }
    }

//    cells[0]._grad_rho = -cells[1]._grad_rho;
//    cells[0]._grad_u = -cells[1]._grad_u - 4. * cells[0]._u / (cells[0]._midpoint - cells[1]._midpoint);
//    cells[0]._grad_P = -cells[1]._grad_P;
    cells[0]._grad_rho = cells[1]._grad_rho;
    cells[0]._grad_u = cells[1]._grad_u;
    cells[0]._grad_P = cells[1]._grad_P;
#if BOUNDARIES == BOUNDARIES_OPEN || BOUNDARIES == BOUNDARIES_BONDI
    cells[NCELL+1]._grad_rho = cells[NCELL+1]._grad_rho;
    cells[NCELL+1]._grad_u = cells[NCELL+1]._grad_u;
    cells[NCELL+1]._grad_P = cells[NCELL+1]._grad_P;
#elif BOUNDARIES == BOUNDARIES_REFLECTIVE
    cells[NCELL+1]._grad_rho = -cells[NCELL]._grad_rho;
    cells[NCELL+1]._grad_u = -cells[NCELL]._grad_u - 4. * cells[NCELL]._u /
                               (cells[NCELL]._midpoint - cells[NCELL+1]._midpoint);
    cells[NCELL+1]._grad_P = -cells[NCELL]._grad_P;
#else
#error "No boundary conditions chosen!"
#endif

#ifdef NO_GRADIENTS
    for(unsigned int i = 0; i < NCELL+2; ++i){
      cells[i]._grad_rho = 0.;
      cells[i]._grad_u = 0.;
      cells[i]._grad_P = 0.;
    }
#endif

    for(unsigned int i = 0; i < NCELL+1; ++i){
      double rhoL, uL, PL, rhoR, uR, PR;
      rhoL = cells[i]._rho;
      uL = cells[i]._u;
      PL = cells[i]._P;
      Cell *rcell = cells[i]._right_ngb;
      rhoR = rcell->_rho;
      uR = rcell->_u;
      PR = rcell->_P;

      double dmin, dplu;
      dmin = 0.5 * (rcell->_midpoint - cells[i]._midpoint);
      dplu = -dmin;
      double rhoL_dash, rhoR_dash, uL_dash, uR_dash, PL_dash, PR_dash;
      rhoL_dash = rhoL + dmin * cells[i]._grad_rho;
      uL_dash = uL + dmin * cells[i]._grad_u;
      PL_dash = PL + dmin * cells[i]._grad_P;
      rhoR_dash = rhoR + dplu * rcell->_grad_rho;
      uR_dash = uR + dplu * rcell->_grad_u;
      PR_dash = PR + dplu * rcell->_grad_P;

      rhoL_dash -= 0.5 * DT * (rhoL * cells[i]._grad_u +
                               uL * cells[i]._grad_rho);
      uL_dash -= 0.5 * DT * (uL * cells[i]._grad_u + cells[i]._grad_P / rhoL);
      PL_dash -= 0.5 * DT * (GAMMA * PL * cells[i]._grad_u +
                             uL * cells[i]._grad_P);
      rhoR_dash -= 0.5 * DT * (rhoR * rcell->_grad_u + uR * rcell->_grad_rho);
      uR_dash -= 0.5 * DT * (uR * rcell->_grad_u + rcell->_grad_P / rhoR);
      PR_dash -= 0.5 * DT * (GAMMA * PR * rcell->_grad_u + uR * rcell->_grad_P);

#if POTENTIAL != POTENTIAL_NONE
      uL_dash += 0.5 * DT * cells[i]._a;
      uR_dash += 0.5 * DT * rcell->_a;
#endif

      double rhosol, usol, Psol;
      solver.solve(rhoL_dash, uL_dash, PL_dash, rhoR_dash, uR_dash, PR_dash,
                   rhosol, usol, Psol);
      double mflux, pflux, Eflux;
      mflux = rhosol * usol;
      pflux = rhosol * usol * usol + Psol;
      Eflux = (Psol / (GAMMA - 1.) + 0.5 * rhosol * usol * usol) * usol +
              Psol * usol;
      cells[i]._m -= DT * mflux;
      cells[i]._p -= DT * pflux;
      cells[i]._E -= DT * Eflux;
      rcell->_m += DT * mflux;
      rcell->_p += DT * pflux;
      rcell->_E += DT * Eflux;
    }

    // add spherical source term
    for(unsigned int i = 1; i < NCELL+1; ++i){
      const double r = cells[i]._midpoint;
      const double Ui[3] = {cells[i]._m / cells[i]._V,
                            cells[i]._p / cells[i]._V,
                            cells[i]._E / cells[i]._V};
      double K1[3];
      K1[0] = -DT * Ui[1] / r;
      K1[1] = -DT * Ui[1] * Ui[1] / Ui[0] / r;
      const double p1 = (GAMMA - 1.) * (Ui[2] - 0.5 * Ui[1] * Ui[1] / Ui[0]);
      K1[2] = -DT * Ui[1] * (Ui[2] + p1) / Ui[0] / r;
      double U[3] = {Ui[0] + K1[0],
                     Ui[1] + K1[1],
                     Ui[2] + K1[2]};
      double K2[3];
      K2[0] = -DT * U[1] / r;
      K2[1] = -DT * U[1] * U[1] / U[0] / r;
      const double p2 = (GAMMA - 1.) * (U[2] - 0.5 * U[1] * U[1] / U[0]);
      K2[2] = -DT * U[1] * (U[2] + p2) / U[0] / r;
      U[0] = Ui[0] + 0.5 * (K1[0] + K2[0]);
      U[1] = Ui[1] + 0.5 * (K1[1] + K2[1]);
      U[2] = Ui[2] + 0.5 * (K1[2] + K2[2]);

      cells[i]._m = U[0] * cells[i]._V;
      cells[i]._p = U[1] * cells[i]._V;
      cells[i]._E = U[2] * cells[i]._V;
    }

    // add gravitational acceleration
#if POTENTIAL == POTENTIAL_POINT_MASS
    for(unsigned int i = 1; i < NCELL+1; ++i){
      const double r = cells[i]._midpoint;
      const double a = -G * MASS_POINT_MASS / (r * r);
      cells[i]._a = a;
      const double m = cells[i]._V * cells[i]._rho;
      cells[i]._p += 0.5 * DT * a * m;
    }
#endif

  }

  write_snapshot(NSTEP, cells);
  for(unsigned int i = 1; i < NCELL+1; ++i){
    std::cout << cells[i]._midpoint << "\t" << cells[i]._rho << "\t"
              << cells[i]._u << "\t" << cells[i]._P << std::endl;
  }

  return 0;
}
