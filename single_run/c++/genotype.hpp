#ifndef GENOTYPE
#define GENOTYPE

#include <vector>
#include <cmath>
#include <iostream>

class Genotype{
private:
  // type-G inheritated from a male-like parent
  std::vector<double> genotype1;
  // type-G inheritated from a female-like parent
  std::vector<double> genotype2;
  // type-B inheritated from a male-like parent
  std::vector<std::vector<double>> interactions1;
  // type-G inheritated from a female-like parent
  std::vector<std::vector<double>> interactions2;
  int num_g;

  // the number of individuals who have this genotype
  int ind_num;
  // phenotypes during the adult phase
  std::vector<std::vector<double>> phenotypes;

public:
  Genotype(const int input_ind_num,
    const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2,
    const std::vector<std::vector<double>>& input_interactions1,
    const std::vector<std::vector<double>>& input_interactions2);
  void set_genotype(const int input_ind_num,
    const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2,
    const std::vector<std::vector<double>>& input_interactions1,
    const std::vector<std::vector<double>>& input_interactions2);
  void development(const double alpha, const int round, const int adult_step);
  std::vector<double> ret_fitness(const double alpha, const int round,
    const int adult_period, const std::vector<std::vector<double>>& s_vecs, 
    const double deno);
  void add_mut(const int index, const double delta, 
    const double max_g, const bool parent);
  void return_genotype1(std::vector<double>& ret) const{ ret = genotype1; };
  void return_genotype2(std::vector<double>& ret) const{ ret = genotype2; };
  void return_interactions1(std::vector<std::vector<double>>& ret) const
    { ret = interactions1; };
  void return_interactions2(std::vector<std::vector<double>>& ret) const
    { ret = interactions2; };
  void return_phenotype(std::vector<double>& ret) const{ ret = phenotypes.back(); };
  bool check_identical(const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2,
    const std::vector<std::vector<double>>& input_interactions1,
    const std::vector<std::vector<double>>& input_interactions2) const;
  void add_ind(){ ind_num++; };
  void reduce_ind(){ ind_num--; };
  void set_ind_num(const int num){ ind_num = num; };
  int ret_ind_num() const{ return(ind_num); };
};

#endif