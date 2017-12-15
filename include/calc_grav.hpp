#pragma once

void clear_acc_and_phi(real_particle * bp, int n);
void accumulate_force_from_point(vector dx, real r2, real eps2,
                        				 vector & acc, real & phi,real jmass);
void tree_walk_counter(real_particle * rp, bhparticle * bp, int nbody, vector cmpos, real length, int level);
void clear_tree_counters();
void print_tree_counters();
void calculate_gravity_using_tree(real_particle * rp, bhnode * bn, real eps2, real theta2, int n_th);
void calculate_gravity(int nbody, real_system pb, bhnode * bn, bhparticle * bp, real eps2, real theta2);
void accumulate_mutual_gravity(real_particle & p1, real_particle & p2, real eps2);
void calculate_gravity_direct(real_particle * pb, int n, double eps2);
