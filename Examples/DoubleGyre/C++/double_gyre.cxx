
#include <math.h>
#include <iostream>
#include <string.h>
#include <string>
#include <algorithm>
#include <cassert>
#include <numbers>

#include "double_gyre.h"

namespace double_gyre
{

void
DoubleGyre::AllocateGrid(int xresolution, int yresolution)
{
  this->timestep = 0.1;

  this->xres = xresolution; // X horizontal resolution
  this->yres = yresolution; // Y vertical   resolution
  assert (this->xres == 2 * this->yres);
  this->xaxis.resize(this->xres);
        
  float spacing = grid_bounds[1] / (this->xres - 1.0);
  for(auto i=0; i< this->xres; i++){
    this->xaxis[i] = i * spacing;
  }
        
  spacing = grid_bounds[3] / (this->yres - 1.0);
  this->yaxis.resize(this->yres);
  for(auto i=0; i< this->yres; i++){
    this->yaxis[i] = i * spacing;
  }
  this->vel_x.resize(this->xres * this->yres);
  this->vel_y.resize(this->xres * this->yres);
  this->vel_z.resize(this->xres * this->yres);
  std::fill(this->vel_z.begin(), this->vel_z.end(), 0.0);
        
  this->A = 0.1 * std::numbers::pi;
  this->w = 2.0 * std::numbers::pi/10.;
  this->E = 0.25;
}

void DoubleGyre::compute_step(void)
{
// Computes and updates velocity fields
  double At, Bt, Ft, fft;
  At = this->E * sin(this->w * this->iteration * this->timestep);
  Bt = 1.0 - 2.0 * At;
  for(auto iy=0; iy < this->yres; iy++)
    for(auto ix=0; ix < this->xres; ix++)
      {
      Ft = (At * this->xaxis[ix]*this->xaxis[ix] + Bt * this->xaxis[ix]) * std::numbers::pi;
      fft = 2.0 * At * this->xaxis[ix] + Bt;
      this->vel_x[iy*xres + ix] = -this->A * sin(Ft) * cos(std::numbers::pi*this->yaxis[iy]);
      this->vel_y[iy*xres + ix] =  this->A * cos(Ft) * sin(std::numbers::pi*this->yaxis[iy])*fft;
      }
  this->iteration++;
}

DoubleGyre simulation;
}
