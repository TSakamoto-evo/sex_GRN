#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include "population.hpp"
#include "parameters.hpp"
#include "genotype.hpp"

void prepare_files(const std::string seed_list, const std::string outcome, 
  const std::string freq, int& rep);

void copy_state(Population& source, Population& recipient);

int main(int argc, char* argv[]){
  if(argc != 6) {
    std::cout << "Unexpected arguments" << std::endl;
    return 1;
  }

  unsigned int seed = static_cast<unsigned int>(std::stoul(argv[1]));
  int start_gen = std::stoi(argv[2]);
  int anc_locus = std::stoi(argv[3]);
  int der_locus = std::stoi(argv[4]);
  double effect = std::stod(argv[5]);

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

  int max_time = 5000000;

  int max_test = 1000;
  int rep_each_state = 1000;
  int rep_each_state_short = 10000;
  int max_gen = 10;
  int ini_rep = 1;

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

  std::ofstream ofs0("output_seed.txt");
  ofs0 << seed << std::endl;
  ofs0.close();

  prepare_files("seed_list.txt", "outcome.txt", "freq.txt", ini_rep);

  Population pop(para, s_vecs, seed);
  std::ofstream ofs1("seed_list.txt", std::ios::app);

  std::ofstream ofs2("outcome.txt", std::ios::app);
  std::ofstream ofs3("freq.txt", std::ios::app);
  std::ofstream ofs4("fitness.txt");

  // run simulation until reaching the specified generation
  for(int i = 0; i <= 2 * max_time; i++){
    if(pop.return_gen() == start_gen){
      break;
    }

    pop.one_generation();

    if(i % 100 == 0){
      std::vector<double> regi;
      pop.return_fitness(regi);

      ofs4 << i;
      for(const auto& j: regi){
        ofs4 << "\t" << j;
      }
      ofs4 << std::endl;
    }
  }

  // record the current state
  Population record_state = pop;
  std::random_device gen_seed;
  std::mt19937 mt;

  int current = pop.return_gen();

  // run independent replication
  for(int rep = ini_rep; rep <= max_test; rep++){
    // start from the recorded state
    pop = record_state;

    // assign a new seed
    unsigned int seed = gen_seed();
    ofs1 << rep << "\t" << 0 << "\t" << seed << std::endl;

    mt.seed(seed);
    pop.set_new_rand(mt);

    // 0: anc, 1: der, 2: irregular, 3: other, 4: anc_new, 5: der_new
    std::vector<int> count(6);
    std::vector<double> freqs(max_gen + 1);

    for(int rep2 = 1; rep2 <= rep_each_state; rep2++){
      pop.return_rand(mt);
      copy_state(record_state, pop);
      pop.set_new_rand(mt);

      pop.introduce_mutation(der_locus, effect);
      pop.output_freq_site(der_locus, effect, freqs, 0);

      int locus = -1;
      int now_time = 0;
      for(now_time = 0; now_time <= max_time; now_time++){
        pop.one_generation();

        if(now_time + 1 <= max_gen){
          pop.output_freq_site(der_locus, effect, freqs, now_time + 1);
        }

        locus = pop.polymorphic(anc_locus, der_locus);
        if(locus == -2){
          count.at(2)++;
          break;
        }else if(locus >= 0){
          break;
        }

        if(now_time == max_time){
          count.at(2)++;
          break;
        }
      }

      if(locus >= 0){
        int state = pop.return_state(anc_locus, der_locus, current);
        count.at(state)++;
      }
    }

    for(int rep2 = rep_each_state + 1; rep2 <= rep_each_state_short; rep2++){
      pop.return_rand(mt);
      copy_state(record_state, pop);
      pop.set_new_rand(mt);

      pop.introduce_mutation(der_locus, effect);
      pop.output_freq_site(der_locus, effect, freqs, 0);

      for(int i = 0; i < max_gen; i++){
        pop.one_generation();
        bool state = pop.output_freq_site(der_locus, effect, freqs, i + 1);

        if(state == 0){
          break;
        }
      }
    }

    for(int i = 0; i <= max_gen; i++){
      ofs3 << rep << "\t" << i << "\t" << freqs.at(i) / rep_each_state_short << "\n"; 
    }
    ofs3 << std::flush;

    ofs2 << rep << "\t" << count.at(0) << "\t" << count.at(1) << "\t" << count.at(2)
      << "\t" << count.at(3) << "\t" << count.at(4) << "\t" << count.at(5) << std::endl;
  }

  return(0);
}

void prepare_files(const std::string seed_list, const std::string outcome, 
  const std::string freq, int& rep){

  if(std::filesystem::exists(seed_list) && std::filesystem::exists(outcome) && 
    std::filesystem::exists(freq)){

    int last_rep = -1;

    {
      // read outcome
      std::ifstream ifs(outcome);
      if(!ifs){
        std::cerr << "Fail to open the outcome file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        last_rep = std::stoi(list.at(0));
      }
    }
    rep = last_rep + 1;

    std::ofstream ofs1("tmp_" + outcome);
    // outcome
    {
      std::ifstream ifs(outcome);
      if(!ifs){
        std::cerr << "Fail to open the outcome file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(std::stoi(list.at(0)) <= last_rep){
          ofs1 << line << "\n";
        }
      }
    }

    std::ofstream ofs2("tmp_" + seed_list);
    // seed list
    {
      std::ifstream ifs(seed_list);
      if(!ifs){
        std::cerr << "Fail to open the seedlist file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(std::stoi(list.at(0)) <= last_rep){
          ofs2 << line << "\n";
        }
      }
    }

    std::ofstream ofs3("tmp_" + freq);

    // freq
    {
      std::ifstream ifs(freq);
      if(!ifs){
        std::cerr << "Fail to open the freq file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(std::stoi(list.at(0)) <= last_rep){
          ofs3 << line << "\n";
        }
      }
    }

    ofs1.close();
    ofs2.close();
    ofs3.close();

    std::filesystem::rename("tmp_" + outcome, outcome);
    std::filesystem::rename("tmp_" + seed_list, seed_list);
    std::filesystem::rename("tmp_" + freq, freq);
  }else{
    std::ofstream ofs1(seed_list);
    std::ofstream ofs2(outcome);
    std::ofstream ofs3(freq);

    rep = 1;
  }
}

void copy_state(Population& source, Population& recipient){
  const std::vector<Genotype>& source_pop = source.ref_pop();
  const int gen = source.return_gen();
  const int max_index = source.return_max_index();
  const std::vector<std::unordered_map<double, int>>& source_mut_origin = source.ref_mut_origin();

  std::vector<Genotype>& recipient_pop = recipient.ref_pop();
  std::vector<std::unordered_map<double, int>>& recipient_pop_mut_origin = recipient.ref_mut_origin();

  for(int i = 0; i <= max_index; i++){
    recipient_pop.at(i) = source_pop.at(i);
  }

  recipient_pop_mut_origin = source_mut_origin;
  recipient.set_gen_and_max_index(gen, max_index);
}