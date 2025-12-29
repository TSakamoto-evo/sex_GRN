#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include "population.hpp"
#include "parameters.hpp"

int main(){
  Parameters para;

  /* Parameter specification */
  // diploid population size
  para.pop_size = 1000;
  // number of genes
  para.num_g = 8;
  // decay speed of gene product
  para.alpha = 0.2;

  // maximum effect size of type-G mutations
  para.gamma_g = 0.1;
  // type-G mutation rate (low: 2e-7, default: 1e-6, high: 5e-6)
  para.mu_g = 1e-6;
  // maximum effect size of type-B mutations
  para.gamma_b = 0.1;
  // type-B mutation rate (low: 2e-7, default: 1e-6, high: 5e-6)
  para.mu_b = 1e-6;

  // selection strength (strong: 3.0, default: 5.0, weak: 7.0)
  para.sig = 5.0;

  // criteria for well-evolved sexes
  para.criteria = 0.95;

  // total developmental steps 
  para.dev_period = 40;
  // length of adult steps
  para.adult_step = 20;

  // fitness reduction caused by a introduced deleterious mutation
  para.s = 0.9;

  // optimum value for high expression genes
  double high = 10.0;
  // optimum value for low expression genes
  double low = 0.0;

  // resolution for the record of phenotype distribution
  double bin = 0.01;

  // maximum waiting time for the evolution of sex
  int first_time = 5000000;
  // time for stabilizing a sex-determining system
  int stabilize_time = 500000;
  // maximum length after the introduction of the deleterious mutation
  int max_time = 5000000;

  // optimum expression profile in each sex
  // 0: low in both sexes, 1: high in both sexes
  // 2: high in male and low in female
  // -2: low in male and high in female

  // default setting
  std::vector<int> tmp_s = {1, 1, 0, 0, 2, 2, -2, -2};
  
  // divergent setting
  //std::vector<int> tmp_s = {2, 2, 2, 2, -2, -2, -2, -2};

  /* End parameter specification */

  para.max_g = 1.0 / para.alpha;

  std::vector<double> s1_vec;
  std::vector<double> s2_vec;

  for(const auto& i: tmp_s){
    if(i == 0){
      s1_vec.push_back(low);
      s2_vec.push_back(low);
    }else if(i == 1){
      s1_vec.push_back(high);
      s2_vec.push_back(high);
    }else if(i == 2){
      s1_vec.push_back(high);
      s2_vec.push_back(low);
    }else if(i == -2){
      s1_vec.push_back(low);
      s2_vec.push_back(high);
    }
  }

  std::vector<std::vector<double>> s_vecs = {s1_vec, s2_vec};

  // seed generation
  std::random_device gen_seed;
  unsigned int seed = gen_seed();

  // record a seed value
  std::ofstream ofs0("output_seed.txt");
  ofs0 << seed << std::endl;
  ofs0.close();

  // initialization
  Population pop(para, s_vecs, seed);

  /* Output files */
  // mean fitness (max(Wm, Wf)), mean fitness of male-like and fenale-like parent across mating pairs, 
  // proportion of each sex, the number of survived adults who enters a mating pool,
  // and frequency of the deleterious allele in each sex
  std::ofstream ofs1("fitness.txt");
  // phenotypic distribution
  std::ofstream ofs2("regi_binned.txt");
  // mean and variance of genotypic values at each locus
  std::ofstream ofs3("mean_var.txt");
  // record timings of the key events
  std::ofstream ofs6("regi_time.txt");
  // association between genotypic values and maleness scores
  std::ofstream ofs7("log_association.txt");
  // genotypic value, origin, frequency of each allele
  std::ofstream ofs8("output_freq.txt");

  // check if the two sexes have evolved
  int success = 0;
  int generation = -1;

  // until the establishment of sex
  for(int i = 0; i <= first_time - 1; i++){
    if(i % 10000 == 0){
      std::vector<double> means, vars;
      pop.return_mean_vars(means, vars);

      for(int j = 0; j < static_cast<int>(means.size()); j++){
        ofs3 << i << "\t" << j << "\t" << means.at(j) << "\t" << vars.at(j)
          << std::endl;
      }
    }

    if(i % 10000 == 0){
      std::vector<int> gene_index, bin_index, ret_all_ind, ret_mm, ret_mf, ret_fm, ret_ff;
      pop.return_binned_distribution(bin, gene_index, bin_index, ret_all_ind,
        ret_mm, ret_mf, ret_fm, ret_ff);

      for(int j = 0; j < static_cast<int>(gene_index.size()); j++){
        ofs2 << i << "\t" << gene_index.at(j) << "\t" << bin_index.at(j) * bin << "\t" <<
          (bin_index.at(j) + 1) * bin << "\t" << ret_all_ind.at(j) << "\t" << ret_mm.at(j) << "\t" <<
          ret_mf.at(j) << "\t" << ret_fm.at(j) << "\t" << ret_ff.at(j) << std::endl;
      }
    }

    if(i % 10000 == 0){
      std::vector<double> ret;
      pop.return_log_association(ret);

      ofs7 << i;
      for(const auto& j: ret){
        ofs7 << "\t" << j;
      }
      ofs7 << std::endl;
    }

    pop.one_generation();

    if(i % 1000 == 0){
      std::vector<double> regi;
      pop.return_fitness(regi);

      ofs1 << i;
      for(const auto& j: regi){
        ofs1 << "\t" << j;
      }
      ofs1 << std::endl;
    }

    if(success == 0 && pop.return_dimorphic()){
      ofs6 << i << std::endl;
      pop.output_freq_origin(ofs8);

      success = 1;
      generation = i;
      break;
    }
  }

  if(generation == -1){
    ofs6 << "NA" << std::endl;
    std::exit(1);
  }

  // stabilizing time
  for(int i = generation + 1; i <= generation + stabilize_time; i++){
    if(i % 10000 == 0){
      std::vector<double> means, vars;
      pop.return_mean_vars(means, vars);

      for(int j = 0; j < static_cast<int>(means.size()); j++){
        ofs3 << i << "\t" << j << "\t" << means.at(j) << "\t" << vars.at(j)
          << std::endl;
      }
    }

    if(i % 10000 == 0){
      std::vector<int> gene_index, bin_index, ret_all_ind, ret_mm, ret_mf, ret_fm, ret_ff;
      pop.return_binned_distribution(bin, gene_index, bin_index, ret_all_ind,
        ret_mm, ret_mf, ret_fm, ret_ff);

      for(int j = 0; j < static_cast<int>(gene_index.size()); j++){
        ofs2 << i << "\t" << gene_index.at(j) << "\t" << bin_index.at(j) * bin << "\t" <<
          (bin_index.at(j) + 1) * bin << "\t" << ret_all_ind.at(j) << "\t" << ret_mm.at(j) << "\t" <<
          ret_mf.at(j) << "\t" << ret_fm.at(j) << "\t" << ret_ff.at(j) << std::endl;
      }
    }

    if(i % 10000 == 0){
      std::vector<double> ret;
      pop.return_log_association(ret);

      ofs7 << i;
      for(const auto& j: ret){
        ofs7 << "\t" << j;
      }
      ofs7 << std::endl;
    }

    pop.one_generation();

    if(i % 1000 == 0 || i == generation + stabilize_time){
      std::vector<double> regi;
      pop.return_fitness(regi);

      ofs1 << i;
      for(const auto& j: regi){
        ofs1 << "\t" << j;
      }
      ofs1 << std::endl;
    }
  }

  pop.output_freq_origin(ofs8);

  success = 0;
  // deleterious mutation
  pop.introduce_del();

  int turnover = -1;

  // after the introduction of the deleterious mutation
  for(int i = generation + stabilize_time + 1; i <= generation + stabilize_time + max_time; i++){
    if((i - generation - stabilize_time - 1) % 10000 == 0 || i == turnover + 1){
      std::vector<double> means, vars;
      pop.return_mean_vars(means, vars);

      for(int j = 0; j < static_cast<int>(means.size()); j++){
        ofs3 << i << "\t" << j << "\t" << means.at(j) << "\t" << vars.at(j)
          << std::endl;
      }
    }

    if((i - generation - stabilize_time - 1) % 10000 == 0 || i == turnover + 1){
      std::vector<int> gene_index, bin_index, ret_all_ind, ret_mm, ret_mf, ret_fm, ret_ff;
      pop.return_binned_distribution(bin, gene_index, bin_index, ret_all_ind,
        ret_mm, ret_mf, ret_fm, ret_ff);

      for(int j = 0; j < static_cast<int>(gene_index.size()); j++){
        ofs2 << i << "\t" << gene_index.at(j) << "\t" << bin_index.at(j) * bin << "\t" <<
          (bin_index.at(j) + 1) * bin << "\t" << ret_all_ind.at(j) << "\t" << ret_mm.at(j) << "\t" <<
          ret_mf.at(j) << "\t" << ret_fm.at(j) << "\t" << ret_ff.at(j) << std::endl;
      }
    }

    if(((i - generation - stabilize_time - 1) < 10000 && (i - generation - stabilize_time - 1) % 100 == 0) || 
      (i - generation - stabilize_time - 1) % 10000 == 0 || i == turnover + 1){

      std::vector<double> ret;
      pop.return_log_association(ret);

      ofs7 << i;
      for(const auto& j: ret){
        ofs7 << "\t" << j;
      }
      ofs7 << std::endl;
    }

    pop.one_generation();

    if(((i - generation - stabilize_time - 1) < 100000 && (i - generation - stabilize_time - 1) % 10 == 0) || 
      (i - generation - stabilize_time - 1) % 10000 == 0 || (success == 0 && pop.return_del_allele() == 0)){

      std::vector<double> regi;
      pop.return_fitness(regi);

      ofs1 << i;
      for(const auto& j: regi){
        ofs1 << "\t" << j;
      }
      ofs1 << std::endl;
    }

    if(success == 0 && pop.return_del_allele() == 0){
      ofs6 << i << std::endl;
      ofs6 << pop.return_old_sd_index() << std::endl;
      success = 1;
      turnover = i;

      pop.output_freq_origin(ofs8);
    }

    if(success == 1 && pop.return_dimorphic()){
      ofs6 << i << std::endl;
      success = 2;

      pop.output_freq_origin(ofs8);
    }
  }

  if(success == 0){
    ofs6 << "NA" << std::endl;
    ofs6 << "NA" << std::endl;
    ofs6 << "NA" << std::endl;
  }else if(success == 1){
    ofs6 << "NA" << std::endl;
  }

  return(0);
}
