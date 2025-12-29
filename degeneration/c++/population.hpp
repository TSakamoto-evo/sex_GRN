#ifndef POPULATION
#define POPULATION

#include "genotype.hpp"
#include "parameters.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <boost/random.hpp>

class Population{
private:
  std::vector<Genotype> pop;
  std::vector<Genotype> next_gen;

  Parameters para;
  std::mt19937 mt;
  int gen;

  std::vector<int> ind_nums;
  std::vector<double> fitness_male;
  std::vector<double> fitness_female;

  std::vector<double> fitness_male2;
  std::vector<double> fitness_female2;

  int count_mut_male;
  int count_mut_female;

  std::vector<std::vector<double>> s_vecs;

  std::uniform_real_distribution<> uni;
  std::poisson_distribution<> mut_num;
  std::discrete_distribution<> choose_type;
  std::uniform_int_distribution<> choose_locus;
  std::uniform_int_distribution<> choose_uni_ind;
  std::uniform_real_distribution<> d_gamma_g;
  std::uniform_real_distribution<> d_gamma_b;
  std::bernoulli_distribution p;

  int max_index;
  int mode;

  int old_sd_index;
  bool dimorphic;
  int del_allele;

  std::vector<std::unordered_map<double, int>> mut_origin;

public:
  Population(const Parameters input_para,
    const std::vector<std::vector<double>>& input_s_vecs,
    const unsigned int seed);
  void development(const bool rand);
  void selection();
  void gen_mutation();
  void one_generation();
  int return_max_index() const{ return(max_index); };
  bool return_dimorphic() const{ return(dimorphic); };
  int return_del_allele() const{ return(del_allele); };
  int return_old_sd_index() const{ return(old_sd_index); };
  void return_fitness(std::vector<double>& regi) const;
  void return_genotypes(std::vector<int>& num_dist,
    std::vector<std::vector<double>>& ret1, std::vector<double>& ret2);
  void return_mean_vars(std::vector<double>& means,
    std::vector<double>& vars) const;
  void return_binned_distribution(const double bin,
    std::vector<int>& gene_index,
    std::vector<int>& bin_index, std::vector<int>& ret_all_ind,
    std::vector<int>& ret_mm, std::vector<int>& ret_mf,
    std::vector<int>& ret_fm, std::vector<int>& ret_ff);
  std::vector<double> calculate_p(const double bin);
  void introduce_del();
  void erase_unused_index();

  double calculate_kendall_tau(const std::vector<double>& list1,
    const std::vector<double>& list2) const;
  double add_log_val(const double val1, const double val2) const;
  double subtract_log_val(const double val1, const double val2) const;
  void return_log_association(std::vector<double>& ret);
  void output_freq_origin(std::ofstream& ofs);
};

#endif
