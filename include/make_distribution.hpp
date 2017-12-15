#pragma once

const vector uniform_random_position_in_sphere(double power_index);
real gaussian(void);
const vector gaussian_velocity(real v_rms);

void shift_theta_hill(int nplummer, double shift_theta[], double shift_radius[]);
void shift_pos_cm(real_particle * bp, int nbody);
void shift_pos_cm_radius(real_particle * bp, int I, int nbody, double r_cut);
void shift_pos_theta_radius(real_particle * bp, int I, int nbody,
                            double shift_theta[], double shift_radius[]);

real randunit(int seed);
real randinter(real a, real b);
const vector uniform_position_in_plummer();
const vector uniform_velocity_in_plummer(vector pos);

void NFW_circular_velocity(real_particle * pb, int nbody);
void NFW_acc(real_particle * pb, int nbody);

void Calculate_center_plummers(real_particle * pb, int nbody, double r_cut, vector cmdx);
void rotating_center_of_mass(real_particle * pb, int nplummer, int nbody, int check,
                              double t, double shift_radius[], double shift_theta[],
                              double center_pos[3]);
void set_initial_distribution(real_system pb, int * N, int * n, int * n_all, double * eps2,
                              double shift_theta[3], double shift_radius[3], double center_pos[3]);
