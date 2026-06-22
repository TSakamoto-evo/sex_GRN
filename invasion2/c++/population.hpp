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

  bool dimorphic;
  std::vector<std::unordered_map<double, int>> mut_origin;

public:
  Population(const Parameters input_para,
    const std::vector<std::vector<double>>& input_s_vecs,
    const unsigned int seed);
  void development();
  void selection();
  void gen_mutation();
  void one_generation();
  int return_max_index() const{ return(max_index); };
  bool return_dimorphic() const{ return(dimorphic); };
  void return_fitness(std::vector<double>& regi) const;
  void return_genotypes(std::vector<int>& num_dist,
    std::vector<std::vector<double>>& ret1,
    std::vector<double>& ret2);
  void return_mean_vars(std::vector<double>& means,
    std::vector<double>& vars) const;
  void return_binned_distribution(const double bin,
    std::vector<int>& gene_index,
    std::vector<int>& bin_index, std::vector<int>& ret_all_ind,
    std::vector<int>& ret_mm, std::vector<int>& ret_mf,
    std::vector<int>& ret_fm, std::vector<int>& ret_ff);
  void erase_unused_index();

  double calculate_kendall_tau(const std::vector<double>& list1,
    const std::vector<double>& list2) const;
  double add_log_val(const double val1, const double val2) const;
  double subtract_log_val(const double val1, const double val2) const;
  void return_log_association(std::vector<double>& ret);
  void output_freq_origin(std::ofstream& ofs, const int rep, const int rep2);

  void set_mutation_rate_zero();
  void set_new_rand(const std::mt19937& input_mt);
  void return_rand(std::mt19937& output_mt);
  int single_polymorphic_locus();
  void introduce_mutation(const int locus, const double effect);
  bool output_freq_site(const int der_locus, const double effect, 
    std::vector<double>& freqs, const int time_index);
  int polymorphic(const int anc_locus, const int der_locus);
  int return_state(const int anc_locus, const int der_locus, const int current);

  std::vector<Genotype>& ref_pop(){ return(pop); };
  int return_gen(){ return(gen); };
  std::vector<std::unordered_map<double, int>>& ref_mut_origin(){ return(mut_origin); };

  void set_gen_and_max_index(const int input_gen, const int input_max_index){
    gen = input_gen;
    max_index = input_max_index;
  }
};

#endif
