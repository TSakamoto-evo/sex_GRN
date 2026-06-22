#include "population.hpp"

Population::Population(const Parameters input_para,
  const std::vector<std::vector<double>>& input_s_vecs,
  const unsigned int seed){

  para = input_para;
  s_vecs = input_s_vecs;

  std::vector<double> input_genotype;
  for(int i = 0; i < para.num_g; i++){
    input_genotype.push_back(para.max_g / 2.0);
  }

  std::vector<std::vector<double>> input_interactions;
  for(int i = 0; i < para.num_g; i++){
    std::vector<double> tmp;
    for(int j = 0; j < para.num_g; j++){
      tmp.push_back(0.0);
    }
    input_interactions.push_back(tmp);
  }

  pop.clear();
  pop.emplace_back(para.pop_size, input_genotype, input_genotype,
    input_interactions, input_interactions);
  max_index = 0;
  mode = 1;

  next_gen.clear();
  for(int i = 0; i < para.pop_size; i++){
    next_gen.emplace_back(0, input_genotype, input_genotype,
      input_interactions, input_interactions);
  }

  std::mt19937 tmp_mt(seed);
  mt = tmp_mt;

  gen = 1;

  std::uniform_real_distribution<> tmp_uni(0.0, 1.0);
  std::poisson_distribution<> tmp_mut_num(2 * para.pop_size * (para.num_g * para.mu_g +
    para.num_g * para.num_g * para.mu_b));
  std::vector<double> prob = {para.num_g * para.mu_g,
    para.num_g * para.num_g * para.mu_b};
  std::discrete_distribution<> tmp_choose_type(prob.begin(), prob.end());
  std::uniform_int_distribution<> tmp_choose_locus(0, para.num_g - 1);
  std::uniform_int_distribution<> tmp_choose_uni_ind(0, para.pop_size - 1);
  std::uniform_real_distribution<> tmp_d_gamma_g(-para.gamma_g, para.gamma_g);
  std::uniform_real_distribution<> tmp_d_gamma_b(-para.gamma_b, para.gamma_b);
  std::bernoulli_distribution tmp_p(0.5);

  uni = tmp_uni;
  mut_num = tmp_mut_num;
  choose_type = tmp_choose_type;
  choose_uni_ind = tmp_choose_uni_ind;
  choose_locus = tmp_choose_locus;
  d_gamma_g = tmp_d_gamma_g;
  d_gamma_b = tmp_d_gamma_b;
  p = tmp_p;

  ind_nums.clear();
  fitness_male.clear();
  fitness_female.clear();

  ind_nums.reserve(para.pop_size);
  fitness_male.reserve(para.pop_size);
  fitness_female.reserve(para.pop_size);

  mut_origin.clear();
  mut_origin = std::vector<std::unordered_map<double, int>>(para.num_g * (para.num_g + 1));

  for(int i = 0; i < para.num_g; i++){
    mut_origin.at(i).emplace(para.max_g / 2.0, 0);
  }

  for(int i = para.num_g; i < para.num_g * (para.num_g + 1); i++){
    mut_origin.at(i).emplace(0.0, 0);
  }
}

void Population::development(){
  double deno = 2.0 * para.sig * para.sig;
  dimorphic = 1;

  fitness_male.clear();
  fitness_female.clear();
  ind_nums.clear();

  std::vector<double> ret;

  for(int i = 0; i <= max_index; i++){
    ret = pop.at(i).ret_fitness(para.alpha, para.dev_period,
      para.adult_step, s_vecs, deno);

    fitness_male.push_back(ret.at(0));
    fitness_female.push_back(ret.at(1));

    ind_nums.push_back(pop.at(i).ret_ind_num());

    if(ret.at(0) <= para.criteria && ret.at(1) <= para.criteria){
      dimorphic = 0;
    }
  }
}

void Population::selection(){
  std::discrete_distribution<> choose_ind(ind_nums.begin(), ind_nums.end());

  std::vector<double> geno1_male, geno2_male, geno1_female, geno2_female;
  std::vector<std::vector<double>> inte1_male, inte2_male, inte1_female, inte2_female;

  double max_male = 0.0;
  double max_female = 0.0;

  for(int i = 0; i <= max_index; i++){
    if(max_male < fitness_male.at(i)){
      max_male = fitness_male.at(i);
    }
    if(max_female < fitness_female.at(i)){
      max_female = fitness_female.at(i);
    }
  }

  max_index = -1;

  double max_fit = max_male * max_female;

  for(int i = 0; i < para.pop_size; i++){
    bool success = 0;

    while(success == 0){
      int parent1 = choose_ind(mt);
      int parent2 = choose_ind(mt);

      double fitness_mf = fitness_male.at(parent1) * fitness_female.at(parent2);
      double fitness_fm = fitness_male.at(parent2) * fitness_female.at(parent1);

      bool pair_order;
      double fitness;

      if(fitness_mf > fitness_fm){
        pair_order = 1;
        fitness = fitness_mf;
      }else{
        pair_order = 0;
        fitness = fitness_fm;
      }

      if(uni(mt) * max_fit < fitness){
        success = 1;

        if(pair_order){
          pop.at(parent1).return_genotype1(geno1_male);
          pop.at(parent1).return_genotype2(geno2_male);
          pop.at(parent1).return_interactions1(inte1_male);
          pop.at(parent1).return_interactions2(inte2_male);

          pop.at(parent2).return_genotype1(geno1_female);
          pop.at(parent2).return_genotype2(geno2_female);
          pop.at(parent2).return_interactions1(inte1_female);
          pop.at(parent2).return_interactions2(inte2_female);
        }else{
          pop.at(parent1).return_genotype1(geno1_female);
          pop.at(parent1).return_genotype2(geno2_female);
          pop.at(parent1).return_interactions1(inte1_female);
          pop.at(parent1).return_interactions2(inte2_female);

          pop.at(parent2).return_genotype1(geno1_male);
          pop.at(parent2).return_genotype2(geno2_male);
          pop.at(parent2).return_interactions1(inte1_male);
          pop.at(parent2).return_interactions2(inte2_male);
        }

        for(int j = 0; j < para.num_g; j++){
          if(geno1_male.at(j) != geno2_male.at(j) && p(mt)){
            geno1_male.at(j) = geno2_male.at(j);
          }
          if(geno1_female.at(j) != geno2_female.at(j) && p(mt)){
            geno1_female.at(j) = geno2_female.at(j);
          }

          for(int k = 0; k < para.num_g; k++){
            if(inte1_male.at(j).at(k) != inte2_male.at(j).at(k) && p(mt)){
              inte1_male.at(j).at(k) = inte2_male.at(j).at(k);
            }
            if(inte1_female.at(j).at(k) != inte2_female.at(j).at(k) && p(mt)){
              inte1_female.at(j).at(k) = inte2_female.at(j).at(k);
            }
          }
        }
      }
    }

    if(mode == 1 || gen % 100 == 0){
      int success = -1;
      // check if the identical genotype already exists in the nextgen pop
      for(int j = 0; j <= max_index; j++){
        if(next_gen.at(j).check_identical(geno1_male, geno1_female,
          inte1_male, inte1_female)){
          success = j;
          break;
        }
      }

      if(success >= 0){
        next_gen.at(success).add_ind();
      }else{
        max_index++;
        next_gen.at(max_index).set_genotype(1, geno1_male, geno1_female,
          inte1_male, inte1_female);
      }
    }else{
      max_index++;
      next_gen.at(max_index).set_genotype(1, geno1_male, geno1_female,
        inte1_male, inte1_female);
    }
  }

  if(max_index < 1000){
    mode = 1;
  }else{
    mode = 0;
  }

  for(int i = max_index + 1; i < para.pop_size; i++){
    next_gen.at(i).set_ind_num(0);
  }

  pop = next_gen;
}

void Population::gen_mutation(){
  std::vector<int> ind_num_dist;
  for(int i = 0; i <= max_index; i++){
    ind_num_dist.push_back(pop.at(i).ret_ind_num());
  }

  int num_mut = mut_num(mt);

  for(int i = 0; i < num_mut; i++){
    int ind = choose_uni_ind(mt);

    int index = 0;
    while(ind >= ind_num_dist.at(index)){
      ind -= ind_num_dist.at(index);
      index++;
    }

    if(ind_num_dist.at(index) > 1){
      ind_num_dist.at(index)--;
      ind_num_dist.push_back(1);

      pop.at(index).reduce_ind();
      max_index++;
      pop.at(max_index) = pop.at(index);
      pop.at(max_index).set_ind_num(1);
      index = max_index;
    }

    int type = choose_type(mt);

    if(type == 0){
      int locus = choose_locus(mt);
      double ret = pop.at(index).add_mut(locus, d_gamma_g(mt), para.max_g, p(mt));

      mut_origin.at(locus)[ret] = gen;
    }else if(type == 1){
      int locus1 = choose_locus(mt);
      int locus2 = choose_locus(mt);
      int locus = para.num_g * locus1 + locus2 + para.num_g;

      double ret = pop.at(index).add_mut(locus, d_gamma_b(mt), para.max_g, p(mt));
      mut_origin.at(locus)[ret] = gen;
    }
  }
}

void Population::erase_unused_index(){
  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  std::vector<std::unordered_set<double>> tmp(para.num_g * (para.num_g + 1));

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    for(int j = 0; j < para.num_g; j++){
      if(tmp.at(j).count(tmp_g1.at(j)) == 0){
        tmp.at(j).emplace(tmp_g1.at(j));
      }
      if(tmp.at(j).count(tmp_g2.at(j)) == 0){
        tmp.at(j).emplace(tmp_g2.at(j));
      }
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int locus_index = para.num_g * j + k + para.num_g;

        if(tmp.at(locus_index).count(tmp_i1.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i1.at(j).at(k));
        }
        if(tmp.at(locus_index).count(tmp_i2.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i2.at(j).at(k));
        }
      }
    }
  }

  for(int i = 0; i < para.num_g * (para.num_g + 1); i++){
    std::vector<double> erase_list;
    for(auto& j: mut_origin.at(i)){
      if(tmp.at(i).count(j.first) == 0){
        erase_list.push_back(j.first);
      }
    }

    for(const auto& j: erase_list){
      mut_origin.at(i).erase(j);
    }
  }
}

void Population::one_generation(){
  development();
  selection();
  erase_unused_index();
  gen_mutation();

  gen++;
}

void Population::return_fitness(std::vector<double>& regi) const{
  double fit = 0.0;
  double count1 = 0.0;
  double count2 = 0.0;
  double count3 = 0.0;

  for(int i = 0; i < static_cast<int>(ind_nums.size()); i++){
    if(fitness_male.at(i) > 0.95){
      count1 += ind_nums.at(i);
    }else if(fitness_female.at(i) > 0.95){
      count2 += ind_nums.at(i);
    }else{
      count3 += ind_nums.at(i);
    }

    if(fitness_male.at(i) > fitness_female.at(i)){
      fit += fitness_male.at(i) * ind_nums.at(i);
    }else{
      fit += fitness_female.at(i) * ind_nums.at(i);
    }
  }

  std::vector<double> ret = {fit / para.pop_size, count1 / para.pop_size,
    count2 / para.pop_size, count3 / para.pop_size};

  regi = ret;
}

void Population::return_genotypes(std::vector<int>& num_dist,
  std::vector<std::vector<double>>& ret1, std::vector<double>& ret2){

  num_dist.clear();
  ret1.clear();
  ret2.clear();

  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  development();

  for(int i = 0; i <= max_index; i++){
    std::vector<double> tmp = {fitness_male.at(i), fitness_female.at(i)};

    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    for(const auto& j: tmp_g1){
      tmp.push_back(j);
    }
    for(const auto& j: tmp_g2){
      tmp.push_back(j);
    }
    for(const auto& j: tmp_i1){
      for(const auto& k: j){
        tmp.push_back(k);
      }
    }
    for(const auto& j: tmp_i2){
      for(const auto& k: j){
        tmp.push_back(k);
      }
    }

    ret1.push_back(tmp);
    num_dist.push_back(pop.at(i).ret_ind_num());

    std::vector<double> tmp2;
    pop.at(i).return_phenotype(tmp2);

    double score = 0.0;
    for(int j = 0; j < para.num_g; j++){
      score += (s_vecs.at(0).at(j) - s_vecs.at(1).at(j)) * tmp2.at(j);
    }
    ret2.push_back(score);
  }
}

void Population::return_mean_vars(std::vector<double>& means,
  std::vector<double>& vars) const{

  std::vector<double> tmp_sum(para.num_g * (para.num_g + 1), 0);
  std::vector<double> tmp_sq_sum(para.num_g * (para.num_g + 1), 0);

  std::vector<double> g1, g2;
  std::vector<std::vector<double>> i1, i2;

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(g1);
    pop.at(i).return_genotype2(g2);
    pop.at(i).return_interactions1(i1);
    pop.at(i).return_interactions2(i2);

    int num_size = pop.at(i).ret_ind_num();

    for(int j = 0; j < para.num_g; j++){
      tmp_sum.at(j) += (g1.at(j) + g2.at(j)) * num_size;
      tmp_sq_sum.at(j) += (g1.at(j) * g1.at(j) + g2.at(j) * g2.at(j)) * num_size;
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int index = para.num_g * (j + 1) + k;

        tmp_sum.at(index) += (i1.at(j).at(k) + i2.at(j).at(k)) * num_size;
        tmp_sq_sum.at(index) += (i1.at(j).at(k) * i1.at(j).at(k)
                              + i2.at(j).at(k) * i2.at(j).at(k)) * num_size;
      }
    }
  }

  for(int i = 0; i < para.num_g * (para.num_g + 1); i++){
    tmp_sum.at(i) /= (2.0 * para.pop_size);
    tmp_sq_sum.at(i) /= (2.0 * para.pop_size);

    tmp_sq_sum.at(i) -= tmp_sum.at(i) * tmp_sum.at(i);
  }

  means = tmp_sum;
  vars = tmp_sq_sum;
}

void Population::return_binned_distribution(const double bin,
  std::vector<int>& gene_index,
  std::vector<int>& bin_index, std::vector<int>& ret_all_ind,
  std::vector<int>& ret_mm, std::vector<int>& ret_mf,
  std::vector<int>& ret_fm, std::vector<int>& ret_ff){

  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  development();

  std::vector<std::unordered_map<int, int>> all_map(para.num_g * (para.num_g + 1), std::unordered_map<int, int>());
  std::vector<std::unordered_map<int, int>> mm(para.num_g * (para.num_g + 1), std::unordered_map<int, int>());
  std::vector<std::unordered_map<int, int>> mf(para.num_g * (para.num_g + 1), std::unordered_map<int, int>());
  std::vector<std::unordered_map<int, int>> fm(para.num_g * (para.num_g + 1), std::unordered_map<int, int>());
  std::vector<std::unordered_map<int, int>> ff(para.num_g * (para.num_g + 1), std::unordered_map<int, int>());

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    for(int j = 0; j < para.num_g; j++){
      int index = j;
      int bin_num1 = std::floor(tmp_g1.at(j) / bin);
      int bin_num2 = std::floor(tmp_g2.at(j) / bin);

      if(all_map.at(index).count(bin_num1) == 0){
        all_map.at(index).emplace(bin_num1, 0);
        mm.at(index).emplace(bin_num1, 0);
        mf.at(index).emplace(bin_num1, 0);
        fm.at(index).emplace(bin_num1, 0);
        ff.at(index).emplace(bin_num1, 0);
      }

      if(all_map.at(index).count(bin_num2) == 0){
        all_map.at(index).emplace(bin_num2, 0);
        mm.at(index).emplace(bin_num2, 0);
        mf.at(index).emplace(bin_num2, 0);
        fm.at(index).emplace(bin_num2, 0);
        ff.at(index).emplace(bin_num2, 0);
      }

      all_map.at(index).at(bin_num1) += pop.at(i).ret_ind_num();
      all_map.at(index).at(bin_num2) += pop.at(i).ret_ind_num();

      if(fitness_male.at(i) > 0.95){
        mm.at(index).at(bin_num1) += pop.at(i).ret_ind_num();
        mf.at(index).at(bin_num2) += pop.at(i).ret_ind_num();
      }

      if(fitness_female.at(i) > 0.95){
        fm.at(index).at(bin_num1) += pop.at(i).ret_ind_num();
        ff.at(index).at(bin_num2) += pop.at(i).ret_ind_num();
      }
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int index = (j + 1) * para.num_g + k;
        int bin_num1 = std::floor(tmp_i1.at(j).at(k) / bin);
        int bin_num2 = std::floor(tmp_i2.at(j).at(k) / bin);

        if(all_map.at(index).count(bin_num1) == 0){
          all_map.at(index).emplace(bin_num1, 0);
          mm.at(index).emplace(bin_num1, 0);
          mf.at(index).emplace(bin_num1, 0);
          fm.at(index).emplace(bin_num1, 0);
          ff.at(index).emplace(bin_num1, 0);
        }

        if(all_map.at(index).count(bin_num2) == 0){
          all_map.at(index).emplace(bin_num2, 0);
          mm.at(index).emplace(bin_num2, 0);
          mf.at(index).emplace(bin_num2, 0);
          fm.at(index).emplace(bin_num2, 0);
          ff.at(index).emplace(bin_num2, 0);
        }

        all_map.at(index).at(bin_num1) += pop.at(i).ret_ind_num();
        all_map.at(index).at(bin_num2) += pop.at(i).ret_ind_num();

        if(fitness_male.at(i) > 0.95){
          mm.at(index).at(bin_num1) += pop.at(i).ret_ind_num();
          mf.at(index).at(bin_num2) += pop.at(i).ret_ind_num();
        }

        if(fitness_female.at(i) > 0.95){
          fm.at(index).at(bin_num1) += pop.at(i).ret_ind_num();
          ff.at(index).at(bin_num2) += pop.at(i).ret_ind_num();
        }
      }
    }
  }

  gene_index.clear();
  bin_index.clear();
  ret_all_ind.clear();
  ret_mm.clear();
  ret_mf.clear();
  ret_fm.clear();
  ret_ff.clear();

  for(int i = 0; i < (para.num_g + 1) * para.num_g; i++){
    for(const auto& j: all_map.at(i)){
      gene_index.push_back(i);
      bin_index.push_back(j.first);
      ret_all_ind.push_back(j.second);
      ret_mm.push_back(mm.at(i).at(j.first));
      ret_mf.push_back(mf.at(i).at(j.first));
      ret_fm.push_back(fm.at(i).at(j.first));
      ret_ff.push_back(ff.at(i).at(j.first));
    }
  }
}

double Population::calculate_kendall_tau(const std::vector<double>& list1,
  const std::vector<double>& list2) const{
  // calculate tau b
  int n = static_cast<int>(list1.size());
  double log_v0 = std::log(n) + std::log(n - 1) + std::log(2 * n + 5);

  std::unordered_map<double, int> for_list1, for_list2;
  for(int i = 0; i < n; i++){
    if(for_list1.count(list1.at(i)) == 0){
      for_list1.emplace(list1.at(i), 1);
    }else{
      for_list1.at(list1.at(i))++;
    }

    if(for_list2.count(list2.at(i)) == 0){
      for_list2.emplace(list2.at(i), 1);
    }else{
      for_list2.at(list2.at(i))++;
    }
  }

  double log_vt = -std::numeric_limits<double>::infinity();
  double log_vu = -std::numeric_limits<double>::infinity();
  double log_v11 = -std::numeric_limits<double>::infinity();
  double log_v12 = -std::numeric_limits<double>::infinity();
  double log_v21 = -std::numeric_limits<double>::infinity();
  double log_v22 = -std::numeric_limits<double>::infinity();

  for(const auto& i: for_list1){
    if(i.second >= 2){
      log_vt = add_log_val(log_vt, std::log(i.second) + std::log(i.second - 1) + std::log(2 * i.second + 5));
      log_v11 = add_log_val(log_v11, std::log(i.second) + std::log(i.second - 1));

      if(i.second >= 3){
        log_v12 = add_log_val(log_v12, std::log(i.second) + std::log(i.second - 1) + std::log(i.second - 2));
      }
    }
  }

  for(const auto& i: for_list2){
    if(i.second >= 2){
      log_vu = add_log_val(log_vu, std::log(i.second) + std::log(i.second - 1) + std::log(2 * i.second + 5));
      log_v21 = add_log_val(log_v21, std::log(i.second) + std::log(i.second - 1));

      if(i.second >= 3){
        log_v22 = add_log_val(log_v22, std::log(i.second) + std::log(i.second - 1) + std::log(i.second - 2));
      }
    }
  }

  if(for_list1.size() == 1 || for_list2.size() == 1){
    return(0.0);
  }

  double log_v1 = log_v11 + log_v21 - std::log(n) - std::log(n - 1) - std::log(2.0);
  double log_v2 = log_v12 + log_v22 - std::log(n) - std::log(n - 1) - std::log(n - 2) - std::log(9.0);

  long long int nc = 0;
  long long int nd = 0;

  for(int i = 0; i < n; i++){
    for(int j = i + 1; j < n; j++){
      if(list1.at(i) > list1.at(j) && list2.at(i) > list2.at(j)){
        nc++;
      }else if(list1.at(i) < list1.at(j) && list2.at(i) < list2.at(j)){
        nc++;
      }else if(list1.at(i) > list1.at(j) && list2.at(i) < list2.at(j)){
        nd++;
      }else if(list1.at(i) < list1.at(j) && list2.at(i) > list2.at(j)){
        nd++;
      }
    }
  }

  double log_v = subtract_log_val(add_log_val(log_v0 - std::log(18.0), add_log_val(log_v1, log_v2)),
    add_log_val(log_vt, log_vu) - std::log(18.0));

  double log_zb;
  if(nc > nd){
    log_zb = std::log(nc - nd) - 0.5 * log_v;
  }else{
    log_zb = std::log(nd - nc) - 0.5 * log_v;
  }

  double pi = std::acos(-1.0);

  double log_p_val1 = std::log(0.5 * std::erfc(std::exp(log_zb) / std::sqrt(2.0)));
  double log_p_val2 = -std::exp(2.0 * log_zb) / 2.0 +
    std::log(std::sqrt(2.0 / pi) / std::exp(log_zb) -
    std::sqrt(2.0 / pi) / std::exp(3.0 * log_zb) +
    3.0 * std::sqrt(2.0 / pi) / std::exp(5.0 * log_zb)) - std::log(2.0);

  if(log_zb > 2.5){
    return(log_p_val2 + std::log(2.0));
  }else{
    return(log_p_val1 + std::log(2.0));
  }
}

double Population::add_log_val(const double val1, const double val2) const{
  double ret;

  if(std::isinf(val1) && std::isinf(val2)){
    if(val1 > 0 || val2 > 0){
      ret = std::numeric_limits<double>::infinity();
    }else{
      ret = -std::numeric_limits<double>::infinity();
    }
  }else{
    if(val1 > val2){
      ret = val1 + std::log(1.0 + std::exp(val2 - val1));
    }else{
      ret = val2 + std::log(1.0 + std::exp(val1 - val2));
    }
  }

  return(ret);
}

double Population::subtract_log_val(const double val1, const double val2) const{
  double ret = val1 + std::log(1.0 - std::exp(val2 - val1));
  return(ret);
}

void Population::return_log_association(std::vector<double>& ret){
  ret.clear();
  std::vector<std::vector<double>> tmp_genotype(para.num_g * (para.num_g + 1),
    std::vector<double>(0));
  std::vector<double> scores;

  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  development();

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    std::vector<double> tmp2;
    pop.at(i).return_phenotype(tmp2);

    double score = 0.0;
    for(int j = 0; j < para.num_g; j++){
      score += (s_vecs.at(0).at(j) - s_vecs.at(1).at(j)) * tmp2.at(j);
    }

    int ind_num = pop.at(i).ret_ind_num();

    for(int l = 0; l < ind_num; l++){
      for(int j = 0; j < para.num_g; j++){
        tmp_genotype.at(j).push_back(tmp_g1.at(j) + tmp_g2.at(j));
      }

      for(int j = 0; j < para.num_g; j++){
        for(int k = 0; k < para.num_g; k++){
          tmp_genotype.at((j + 1) * para.num_g + k).push_back(tmp_i1.at(j).at(k) + tmp_i2.at(j).at(k));
        }
      }

      scores.push_back(score);
    }
  }

  for(int i = 0; i < para.num_g * (para.num_g + 1); i++){
    double log_p = calculate_kendall_tau(tmp_genotype.at(i), scores);
    ret.push_back(log_p);
  }
}

void Population::output_freq_origin(std::ofstream& ofs, const int rep, const int rep2){
  std::vector<double> log_p_val;
  return_log_association(log_p_val);

  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  std::vector<std::unordered_map<double, int>> tmp(para.num_g * (para.num_g + 1));

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    int ind_num = pop.at(i).ret_ind_num();

    for(int j = 0; j < para.num_g; j++){
      if(tmp.at(j).count(tmp_g1.at(j)) == 0){
        tmp.at(j).emplace(tmp_g1.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g1.at(j)) += ind_num;
      }
      if(tmp.at(j).count(tmp_g2.at(j)) == 0){
        tmp.at(j).emplace(tmp_g2.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g2.at(j)) += ind_num;
      }
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int locus_index = para.num_g * j + k + para.num_g;

        if(tmp.at(locus_index).count(tmp_i1.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i1.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i1.at(j).at(k)) += ind_num;
        }
        if(tmp.at(locus_index).count(tmp_i2.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i2.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i2.at(j).at(k)) += ind_num;
        }
      }
    }
  }

  for(int i = 0; i < para.num_g * (para.num_g + 1); i++){
    for(const auto j: tmp.at(i)){
      ofs << rep << "\t" << rep2 << "\t" << gen << "\t" << i << "\t" << log_p_val.at(i) << "\t" << j.first << "\t" 
        << mut_origin.at(i).at(j.first) << "\t" << 1.0 * j.second / 2.0 / para.pop_size << "\n";
    }
  }
}

void Population::set_mutation_rate_zero(){
  std::poisson_distribution<> tmp_mut_num(0.0);
  mut_num = tmp_mut_num;

  para.mu_b = 0.0;
  para.mu_g = 0.0;
}

void Population::set_new_rand(const std::mt19937& input_mt){
  mt = input_mt;
}

void Population::return_rand(std::mt19937& output_mt){
  output_mt = mt;
}

int Population::single_polymorphic_locus(){
  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  std::vector<std::unordered_map<double, int>> tmp(para.num_g * (para.num_g + 1));

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    int ind_num = pop.at(i).ret_ind_num();

    for(int j = 0; j < para.num_g; j++){
      if(tmp.at(j).count(tmp_g1.at(j)) == 0){
        tmp.at(j).emplace(tmp_g1.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g1.at(j)) += ind_num;
      }
      if(tmp.at(j).count(tmp_g2.at(j)) == 0){
        tmp.at(j).emplace(tmp_g2.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g2.at(j)) += ind_num;
      }
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int locus_index = para.num_g * j + k + para.num_g;

        if(tmp.at(locus_index).count(tmp_i1.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i1.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i1.at(j).at(k)) += ind_num;
        }
        if(tmp.at(locus_index).count(tmp_i2.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i2.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i2.at(j).at(k)) += ind_num;
        }
      }
    }
  }

  int polymorphic_locus = 0;
  int locus_pos = -1;

  for(int i = 0; i < para.num_g * (para.num_g + 1); i++){
    for(const auto j: tmp.at(i)){
      if(j.second < 2 * para.pop_size){
        polymorphic_locus++;
        locus_pos = i;
        break;
      }
    }
  }

  if(polymorphic_locus >= 2){
    return(-1);
  }else if(polymorphic_locus == 1){
    return(locus_pos);
  }else{
    return(-2);
  }
}

void Population::introduce_mutation(const int locus, const double effect){
  std::vector<int> ind_num_dist;
  for(int i = 0; i <= max_index; i++){
    ind_num_dist.push_back(pop.at(i).ret_ind_num());
  }

  int ind = choose_uni_ind(mt);

  int index = 0;
  while(ind >= ind_num_dist.at(index)){
    ind -= ind_num_dist.at(index);
    index++;
  }

  if(ind_num_dist.at(index) > 1){
    ind_num_dist.at(index)--;
    ind_num_dist.push_back(1);

    pop.at(index).reduce_ind();
    max_index++;
    pop.at(max_index) = pop.at(index);
    pop.at(max_index).set_ind_num(1);
    index = max_index;
  }

  pop.at(index).add_mut_specified(locus, effect, p(mt));
  mut_origin.at(locus)[effect] = gen;
}

bool Population::output_freq_site(const int der_locus, const double effect, 
  std::vector<double>& freqs, const int time_index){
  
  int target_count = 0;

  for(int i = 0; i <= max_index; i++){
    pop.at(i).add_ac(der_locus, effect, target_count);
  }
  freqs.at(time_index) += 1.0 * target_count / 2.0 / para.pop_size;

  return(target_count);
}

int Population::polymorphic(const int anc_locus, const int der_locus){
  
  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  std::vector<std::unordered_map<double, int>> tmp(para.num_g * (para.num_g + 1));

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    int ind_num = pop.at(i).ret_ind_num();

    for(int j = 0; j < para.num_g; j++){
      if(tmp.at(j).count(tmp_g1.at(j)) == 0){
        tmp.at(j).emplace(tmp_g1.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g1.at(j)) += ind_num;
      }
      if(tmp.at(j).count(tmp_g2.at(j)) == 0){
        tmp.at(j).emplace(tmp_g2.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g2.at(j)) += ind_num;
      }
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int locus_index = para.num_g * j + k + para.num_g;

        if(tmp.at(locus_index).count(tmp_i1.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i1.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i1.at(j).at(k)) += ind_num;
        }
        if(tmp.at(locus_index).count(tmp_i2.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i2.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i2.at(j).at(k)) += ind_num;
        }
      }
    }
  }

  size_t anc_size = tmp.at(anc_locus).size();
  size_t der_size = tmp.at(der_locus).size();

  if(anc_size > 1 && der_size == 1){
    return(anc_locus);
  }else if(anc_size == 1 && der_size > 1){
    return(der_locus);
  }else if(anc_size > 1 && der_size > 1){
    return(-1);
  }else{
    return(-2);
  }
}

int Population::return_state(const int anc_locus, const int der_locus, 
  const int current){

  std::vector<double> log_p_val;
  return_log_association(log_p_val);

  int min_index = 0;
  double min_val = log_p_val.at(0);

  for(size_t i = 0; i < log_p_val.size(); i++){
    if(log_p_val.at(i) < min_val){
      min_val = log_p_val.at(i);
      min_index = i;
    }
  }


  std::vector<double> tmp_g1, tmp_g2;
  std::vector<std::vector<double>> tmp_i1, tmp_i2;

  std::vector<std::unordered_map<double, int>> tmp(para.num_g * (para.num_g + 1));

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype1(tmp_g1);
    pop.at(i).return_genotype2(tmp_g2);
    pop.at(i).return_interactions1(tmp_i1);
    pop.at(i).return_interactions2(tmp_i2);

    int ind_num = pop.at(i).ret_ind_num();

    for(int j = 0; j < para.num_g; j++){
      if(tmp.at(j).count(tmp_g1.at(j)) == 0){
        tmp.at(j).emplace(tmp_g1.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g1.at(j)) += ind_num;
      }
      if(tmp.at(j).count(tmp_g2.at(j)) == 0){
        tmp.at(j).emplace(tmp_g2.at(j), ind_num);
      }else{
        tmp.at(j).at(tmp_g2.at(j)) += ind_num;
      }
    }

    for(int j = 0; j < para.num_g; j++){
      for(int k = 0; k < para.num_g; k++){
        int locus_index = para.num_g * j + k + para.num_g;

        if(tmp.at(locus_index).count(tmp_i1.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i1.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i1.at(j).at(k)) += ind_num;
        }
        if(tmp.at(locus_index).count(tmp_i2.at(j).at(k)) == 0){
          tmp.at(locus_index).emplace(tmp_i2.at(j).at(k), ind_num);
        }else{
          tmp.at(locus_index).at(tmp_i2.at(j).at(k)) += ind_num;
        }
      }
    }
  }

  if(min_index == anc_locus){
    int origin1 = -1;
    int origin2 = -1;
    int hapnum1 = 0;
    int hapnum2 = 0;

    for(const auto& i: tmp.at(anc_locus)){
      if(i.second > hapnum1 || (i.second == hapnum1 && origin1 > mut_origin.at(anc_locus).at(i.first))){
        origin2 = origin1;
        hapnum2 = hapnum1;

        origin1 = mut_origin.at(anc_locus).at(i.first);
        hapnum1 = i.second;
      }else if(i.second > hapnum2 || (i.second == hapnum2 && origin2 > mut_origin.at(anc_locus).at(i.first))){
        origin2 = mut_origin.at(anc_locus).at(i.first);
        hapnum2 = i.second;
      }
    }

    if(origin1 > current || origin2 > current){
      return(4);
    }else{
      return(0);
    }
  }else if(min_index == der_locus){
    int origin1 = -1;
    int origin2 = -1;
    int hapnum1 = 0;
    int hapnum2 = 0;

    for(const auto& i: tmp.at(der_locus)){
      if(i.second > hapnum1 || (i.second == hapnum1 && origin1 > mut_origin.at(der_locus).at(i.first))){
        origin2 = origin1;
        hapnum2 = hapnum1;

        origin1 = mut_origin.at(der_locus).at(i.first);
        hapnum1 = i.second;
      }else if(i.second > hapnum2 || (i.second == hapnum2 && origin2 > mut_origin.at(der_locus).at(i.first))){
        origin2 = mut_origin.at(der_locus).at(i.first);
        hapnum2 = i.second;
      }
    }

    if(origin1 > current || origin2 > current){
      return(5);
    }else{
      return(1);
    }
  }else{
    return(3);
  }
}