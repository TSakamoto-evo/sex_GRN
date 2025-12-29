#include <iostream>
#include <fstream>
#include <random>
#include "population.hpp"
#include "parameters.hpp"

int main(){
  Parameters para;
  para.pop_size = 1000;

  para.num_g = 8;
  para.alpha = 0.2;

  para.gamma_g = 0.1;
  para.mu_g = 1e-6;
  para.gamma_b = 0.1;
  para.mu_b = 1e-6;

  para.max_g = 1.0 / para.alpha;
  para.sig = 5.0;

  para.criteria = 0.95;
  para.dev_period = 40;
  para.adult_step = 20;

  double high = 10.0;
  double low = 0.0;

  double bin = 0.01;
  int max_time = 5000000;

  std::vector<int> tmp_s = {1, 1, 0, 0, 2, 2, -2, -2};
  //std::vector<int> tmp_s = {2, 2, 2, 2, -2, -2, -2, -2};

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

  unsigned int seed = 2164286046;

  std::ofstream ofs0("output_seed.txt");
  ofs0 << seed << std::endl;
  ofs0.close();


  Population pop(para, s_vecs, seed);

  std::ofstream ofs1("fitness.txt");
  std::ofstream ofs2("regi_binned.txt");
  std::ofstream ofs3("mean_var.txt");
  std::ofstream ofs4("regi_genotype.txt");
  std::ofstream ofs5("axis_val.txt");
  std::ofstream ofs7("log_association.txt");

  bool success = 0;
  int evo_time = 0;

  for(int i = 0; i <= 2 * max_time; i++){
    if(i % 10000 == 0){
      std::vector<double> means, vars;
      pop.return_mean_vars(means, vars);

      for(int j = 0; j < para.num_g * (para.num_g + 1); j++){
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

    if(i % 10 == 0 && i >= 2315000 && i <= 2335000){
      std::vector<int> num_dist;
      std::vector<std::vector<double>> ret1;
      std::vector<double> ret2;

      pop.return_genotypes(num_dist, ret1, ret2);
      for(int j = 0; j < static_cast<int>(num_dist.size()); j++){
        ofs4 << i << "\t" << j << "\t" << num_dist.at(j);
        for(const auto& k: ret1.at(j)){
          ofs4 << "\t" << k;
        }
        ofs4 << std::endl;

        ofs5 << i << "\t" << ret2.at(j) << std::endl;
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

    if(success == 0 && pop.return_dimorphic()){
      std::ofstream ofs6("regi_time.txt");
      ofs6 << i << std::endl;
      ofs6.close();

      success = 1;
      evo_time = i;
    }

    if(i % 100 == 0){
      std::vector<double> regi;
      pop.return_fitness(regi);

      ofs1 << i;
      for(const auto& j: regi){
        ofs1 << "\t" << j;
      }
      ofs1 << std::endl;
    }

    if((i == max_time && success == 0) || (i == evo_time + max_time && success)){
      break;
    }
  }

  return(0);
}
