#pragma once
class complex{
private:
  double x[2];
public:
  complex(){
    x[0] = 0.0;
    x[1] = 0.0;
  }
  void clear(){
    x[0] = 0.0;
    x[1] = 0.0;
  }
  double get_real() {return x[0];}
  double get_imag() {return x[1];}
  void set_real(double new_real) {x[0] = new_real;}
  void set_imag(double new_imag) {x[1] = new_imag;}
  void set_factor(double factor) {
    x[0] *= factor;
    x[1] *= factor;
  }
  // void set_complex (complex new_comp) {
  //   complex temp = new_comp;
  //   x[0].set_real(temp.get_real());
  //   x[1].set_imag(temp.get_imag());
  // }
  void complex_add (double factor, complex b) {
    x[0] += factor * b.get_real();
    x[1] += factor * b.get_imag();
  }
};
