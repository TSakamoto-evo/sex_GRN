#ifndef PARAMETERS
#define PARAMETERS

class Parameters{
public:
  int pop_size;

  int num_g;
  double alpha;

  double gamma_g;
  double mu_g;
  double gamma_b;
  double mu_b;

  double max_g;

  double sig;
  double criteria;

  int dev_period;
  int adult_step;

  double s;
};

#endif
