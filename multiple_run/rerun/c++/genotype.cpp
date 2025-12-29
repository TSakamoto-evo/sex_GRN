#include "genotype.hpp"

Genotype::Genotype(const int input_ind_num,
  const std::vector<double>& input_genotype1,
  const std::vector<double>& input_genotype2,
  const std::vector<std::vector<double>>& input_interactions1,
  const std::vector<std::vector<double>>& input_interactions2){

  ind_num = input_ind_num;

  genotype1 = input_genotype1;
  genotype2 = input_genotype2;
  interactions1 = input_interactions1;
  interactions2 = input_interactions2;

  num_g = static_cast<int>(genotype1.size());
}

void Genotype::set_genotype(const int input_ind_num,
  const std::vector<double>& input_genotype1,
  const std::vector<double>& input_genotype2,
  const std::vector<std::vector<double>>& input_interactions1,
  const std::vector<std::vector<double>>& input_interactions2){

  ind_num = input_ind_num;

  genotype1 = input_genotype1;
  genotype2 = input_genotype2;
  interactions1 = input_interactions1;
  interactions2 = input_interactions2;

  num_g = static_cast<int>(genotype1.size());
}

void Genotype::development(const double alpha, const int round,
  const int adult_step){

  phenotypes.clear();

  std::vector<double> p_vec(num_g);
  for(int i = 0; i < num_g; i++){
    p_vec.at(i) = genotype1.at(i) + genotype2.at(i);
  }

  for(int i = 0; i < round; i++){
    std::vector<double> tmp_p_vec(num_g);
    for(int j = 0; j < num_g; j++){
      double effect1 = 0.0;
      double effect2 = 0.0;
      for(int k = 0; k < num_g; k++){
        effect1 += interactions1.at(j).at(k) * p_vec.at(k);
        effect2 += interactions2.at(j).at(k) * p_vec.at(k);
      }

      tmp_p_vec.at(j) = (1.0 - alpha) * p_vec.at(j) + 
        (0.5 + 0.5 * std::tanh(effect1)) + (0.5 + 0.5 * std::tanh(effect2));
    }
    p_vec = tmp_p_vec;

    if(i >= round - adult_step){
      phenotypes.push_back(p_vec);
    }
  }
}

std::vector<double> Genotype::ret_fitness(const double alpha, const int round,
  const int adult_step, const std::vector<std::vector<double>>& s_vecs, 
  const double deno){

  development(alpha, round, adult_step);

  std::vector<double> ret;
  for(const auto& s_vec: s_vecs){
    double fitness = 1.0;

    for(const auto& phenotype: phenotypes){
      double distance_sq = 0.0;
      for(int i = 0; i < num_g; i++){
        distance_sq += (s_vec.at(i) - phenotype.at(i)) * (s_vec.at(i) - phenotype.at(i));
      }
      fitness *= std::exp(-distance_sq / deno / adult_step);
    }

    ret.push_back(fitness);
  }
  
  return(ret);
}


void Genotype::add_mut(const int index, const double delta, 
  const double max_g, const bool parent){

  if(parent){
    if(index < num_g){
      genotype1.at(index) += delta;
      if(genotype1.at(index) < 0.0){
        genotype1.at(index) = 0.0;
      }
      if(genotype1.at(index) > max_g){
        genotype1.at(index) = max_g;
      }
    }else{
      int locus1 = (index - num_g) / num_g;
      int locus2 = (index - num_g) % num_g;

      interactions1.at(locus1).at(locus2) += delta;
    }
  }else{
    if(index < num_g){
      genotype2.at(index) += delta;
      if(genotype2.at(index) < 0.0){
        genotype2.at(index) = 0.0;
      }
      if(genotype2.at(index) > max_g){
        genotype2.at(index) = max_g;
      }
    }else{
      int locus1 = (index - num_g) / num_g;
      int locus2 = (index - num_g) % num_g;

      interactions2.at(locus1).at(locus2) += delta;
    }
  }
}

bool Genotype::check_identical(const std::vector<double>& input_genotype1,
  const std::vector<double>& input_genotype2,
  const std::vector<std::vector<double>>& input_interactions1,
  const std::vector<std::vector<double>>& input_interactions2) const{

  for(int i = 0; i < num_g; i++){
    if(genotype1.at(i) != input_genotype1.at(i)){
      return(0);
    }
    if(genotype2.at(i) != input_genotype2.at(i)){
      return(0);
    }
  }

  for(int i = 0; i < num_g; i++){
    for(int j = 0; j < num_g; j++){
      if(interactions1.at(i).at(j) != input_interactions1.at(i).at(j)){
        return(0);
      }
      if(interactions2.at(i).at(j) != input_interactions2.at(i).at(j)){
        return(0);
      }
    }
  }

  return(1);
}