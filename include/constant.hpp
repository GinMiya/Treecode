///******  DEFINE QUANTITIES  *********///
///     Mass of Plummer is 10^8 M_sun  ///
///     HMR of Plummer is 1 kpc        ///
///     Gravitational constant = 1     ///
///     Scale density of NFW profile   ///
///     Scale radius of NFW profile    ///
///   The number of particle of 1bin   ///
///   The number of multipole exp.     ///
///************************************///
#define M0	1.0
#define R0	1.0
#define G0	1.0
#define rho_nfw	(2.23/pow(10.0,2.0))
#define rs_nfw 34.6
#define nbin 100
#define n_threads 8
#define p_max 4
#define delta_int (pow(10.0, -4.0))
