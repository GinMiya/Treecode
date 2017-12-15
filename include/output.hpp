#pragma once
void integrate(int nbody, real_system pb, bhnode * bn, bhparticle * bp,
                real eps2, real theta2, real dt, int check);

int compare_double(const void *a, const void *b);

void compare_acc(int nbody, real_particle * pp, real Sacc[]);
void compare_time(int nbody, real_particle * pp, real Stime[]);

void print_quantities(int nplummer, int n, int n_all, double eps2, double theta,
                      double dt, double tend, int iout, double M, double R,
                      double shift_theta[], double shift_radius[]);
void density_particle(real_particle * pp, int nbody, double SR[], double dens[]);
void search_tree_walk(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2);

void input_number_of_particles(int * n_p, int * N_pl, int * n_all, real shift_radius[]);

void include_quantities(char *arg[], int *N, int *n, int *n_all, double *eps2, double *theta2,
                        double *dt, double *tend,int *iout, double *M, double *R);
void print_simulation_mode(int *N, int *n, int *n_all, double shift_theta[3], double shift_radius[3]);

stringstream make_directory(double shift_radius[3]);

void set_initial_distribution(real_particle * pb, int N, int n, int n_all, double eps2,
                              double shift_theta[3], double shift_radius[3], double center_pos[3]);

void output_acc_error(const char * dirname, real_system pb, int n_all, double Sacc[], double Stime[]);
void output_x_v(const char * dirname, real_system pb, int n_all, int id);
void output_acc(const char * dirname, real_system pb, int n_all, int id);
void output_energy(const char * dirname, real_system pb, int n_all, double t, double center_pos[3]);
void output_dens(const char * dirname, real_system pb, int n_all, int id, double SR[], double dens[]);

void gnuplot_plot();
