%module pyconf
%include <std_vector.i>
%{
  extern void outputPoscar();
  extern void setEquilibriate(double temperature, int num_steps);
  extern const double getTotalEnergy();
  extern const std::vector<int> getSpins();
  extern double metropolis_step(double temperature);
%}


namespace std
{
  %template(IntVector) vector<int>;
}

/* Let's just grab the original header file here */
extern void outputPoscar();
extern void setEquilibriate(double temperature, int num_steps);
extern const double getTotalEnergy();
extern const std::vector<int> getSpins();
extern double metropolis_step(double temperature);
