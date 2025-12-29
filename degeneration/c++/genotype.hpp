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
  // deleterious alleles
  bool allele1;
  bool allele2;

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
    const std::vector<std::vector<double>>& input_interactions2, 
    const bool input_allele1, const bool input_allele2);
  void set_genotype(const int input_ind_num,
    const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2,
    const std::vector<std::vector<double>>& input_interactions1,
    const std::vector<std::vector<double>>& input_interactions2, 
    const bool input_allele1, const bool input_allele2);
  void development(const double alpha, const int round, const int adult_step);
  std::vector<double> ret_fitness(const double alpha, const int round,
    const int adult_period, const std::vector<std::vector<double>>& s_vecs, 
    const double deno);
  double add_mut(const int index, const double delta, 
    const double max_g, const bool parent);
  void introduce_del(const bool a1, const bool a2){ allele1 = a1; allele2 = a2; };
  void return_genotype1(std::vector<double>& ret) const{ ret = genotype1; };
  void return_genotype2(std::vector<double>& ret) const{ ret = genotype2; };
  void return_interactions1(std::vector<std::vector<double>>& ret) const
    { ret = interactions1; };
  void return_interactions2(std::vector<std::vector<double>>& ret) const
    { ret = interactions2; };
  void return_phenotype(std::vector<double>& ret) const{ ret = phenotypes.back(); };
  bool return_allele1() const{ return(allele1); };
  bool return_allele2() const{ return(allele2); };
  bool check_identical(const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2,
    const std::vector<std::vector<double>>& input_interactions1,
    const std::vector<std::vector<double>>& input_interactions2, 
    const bool input_allele1, const bool input_allele2) const;
  void add_ind(){ ind_num++; };
  void reduce_ind(){ ind_num--; };
  void set_ind_num(const int num){ ind_num = num; };
  int ret_ind_num() const{ return(ind_num); };
  int ret_del_num() const{ return(allele1 + allele2); };

  void return_locus(const int index, double& ret1, double& ret2, int& ret_ind_num) const;
  void add_del(const int index, const double healthy);
};

#endif