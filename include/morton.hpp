#pragma once

double ConvertEvenValue(double val);

string to_binString(unsigned int val);

BHlong conv_to_morton_old(int ix, int keybits);
BHlong conv_to_morton(int ix, int keybits);

BHlong dilate3D(BHlong val);
BHlong getMorton3D(const BHlong ix, const BHlong iy, const BHlong iz);
inline BHlong construct_key(const vector & pos, real rscale, int ioffset, int keybits);

int compare_key(const void * p1, const void * p2);

real initialize_key(int nbody, real_particle * rp, int & nkeysize, bhparticle * bhp);
