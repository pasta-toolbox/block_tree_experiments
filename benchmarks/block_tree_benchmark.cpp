/*******************************************************************************
 * This file is part of pasta::block_tree_experiments
 *
 * Copyright (C) 2022-2023 Daniel Meyer
 * Copyright (C) 2023 Florian Kurpicz <florian@kurpicz.org>
 *
 * pasta::block_tree_experiments is free software: you can
 * redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * pasta::block_tree_experiments is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pasta::block_tree_experiments.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <iostream>
#include "libsais.h"
#include <pasta/block_tree/utils/lpf_array.hpp>
#include <vector>
#include <tlx/cmdline_parser.hpp>
#include <fstream>
#include <sstream>
#include <pasta/block_tree/construction/block_tree_lpf.hpp>
#include <pasta/block_tree/construction/block_tree_fp.hpp>
#include <compressed/PCBlockTree.h>
#include <chrono>
#include <type_traits>
#include <iostream>
#include <unordered_set>

class BlockTreeBenchmark {

private:  
  std::string input_path_;
  size_t prefix_size_;
  size_t queries_;

  size_t tau_;
  size_t s_;
  size_t max_leaf_size_;
  
  bool cut_first_levels_;
  bool extended_pruning_;

  size_t threads_;
  bool original_;
  bool parallel_;

  std::vector<uint8_t> input_;
  bool number_lz_phrases_computed_ = false;
  size_t number_lz_phrases_;

  std::vector<uint32_t> access_queries_;
  std::vector<uint32_t> select_queries_;
  std::vector<uint8_t> select_c_;
  
public:

  BlockTreeBenchmark(std::string const& input_path, size_t const prefix_size,
                     size_t const queries, size_t const tau, size_t const s,
                     size_t const max_leaf_size,
                     bool const cut_first_levels, bool const extended_pruning, size_t const threads,
		     bool const original, bool const parallel)
    : input_path_(input_path), prefix_size_(prefix_size), queries_(queries),
      tau_(tau), s_(s), max_leaf_size_(max_leaf_size),
      cut_first_levels_(cut_first_levels), extended_pruning_(extended_pruning), threads_(threads) {

    access_queries_.reserve(queries_);
    select_queries_.reserve(queries_);
    select_c_.reserve(queries_);
   }

  void run() {
    load_input();

    auto const t_lz_start = std::chrono::high_resolution_clock::now();
    compute_number_lz_phrases();
    auto const t_lz_end = std::chrono::high_resolution_clock::now();

    auto const lz_and_z_time = std::chrono::duration_cast<std::chrono::milliseconds>(t_lz_end - t_lz_start).count();
    
    std::cout << "number_lz_phrases() " << number_lz_phrases_ << '\n';

    if (s_ == 0) {
      std::cout << "Tuning parameter s=0: Choosing s=z with " << "\n";
      s_ = number_lz_phrases_;
    }

    prepare_queries();

    if (parallel_) {
      run_construction_lpf(lz_and_z_time, false, threads_);
    } else {
      if (original_) {
	run_construction_original();
      } else {
	run_construction_lpf(lz_and_z_time, true);
	run_construction_lpf(lz_and_z_time, false);
	run_construction_fp(lz_and_z_time);
      }
    }
  }
  
private:

void run_construction_fp(size_t const lz_and_z_time) {

  std::string simple_algo_name = "FP-pruning-" + std::string(s_ == 1 ? "1" : "z");
    
    std::cout << "RESULT"
              << " algorithm=" << simple_algo_name
              << " input=" << input_path_
              << " prefix_size=" << input_.size()
              << " queries=" << queries_
              << " lz_and_z_time=" << lz_and_z_time
              << " s=" << s_
              << " tau=" << tau_
              << " max_leaf_size=" << max_leaf_size_
              << " cut_first_levels=" << (cut_first_levels_ ? "true" : "false")
              << " extended_pruning=" << (extended_pruning_ ? "true" : "false")
              << " dynamic_prog=" << "false";


    auto const t1 = std::chrono::high_resolution_clock::now();
    auto* fpPruned_bt =
      new pasta::BlockTreeFP<uint8_t, int32_t>(input_, tau_, max_leaf_size_,
                                                   s_,
                                                   256 /*sigma, TODO use limit*/,
                                                   cut_first_levels_,
                                                   extended_pruning_);

    auto const space_usage = fpPruned_bt->print_space_usage();

    auto const t2 = std::chrono::high_resolution_clock::now();

    fpPruned_bt->add_rank_support();
    auto const space_usage_rs = fpPruned_bt->print_space_usage();
    
    auto const t3 = std::chrono::high_resolution_clock::now();

    auto const construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    auto const rank_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    std::cout << " space_usage=" << space_usage;
    std::cout << " construction_time_ms=" << construction_time;
    
    std::cout << " space_usage_with_rs=" <<  space_usage_rs;
    std::cout << " construction_time_with_rs_ms=" << construction_time + rank_construction_time;

    size_t result = 0;
    auto const t4 = std::chrono::high_resolution_clock::now();
    for (auto const& query : access_queries_) {
        result += fpPruned_bt->access(query);
    }
    auto const t5 = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < select_c_.size(); i++) {
        result += fpPruned_bt->rank(select_c_[i], access_queries_[i]);
    }

    auto const t6 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < select_c_.size(); i++) {
        result += fpPruned_bt->select(select_c_[i], select_queries_[i]);
    }

    auto const t7 = std::chrono::high_resolution_clock::now();


    auto const access_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count();
    auto const rank_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count();
    auto const select_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t6).count();
    
    std::cout << " access_time_ns=" << access_time
              << " rank_time_ns=" << rank_time
              << " select_time_ns=" << select_time
              << " check_sum=" << result << '\n';
}

  
  void run_construction_lpf(size_t const lz_and_z_time, bool const dynamic_prog, size_t threads = 0) {

    std::string simple_algo_name = "LPF-pruning-" + std::string(dynamic_prog ? "DP-" : "") + (s_ == 1 ? "1" : "z");
    
    std::cout << "RESULT"
              << " algorithm=" << simple_algo_name
              << " input=" << input_path_
              << " prefix_size=" << input_.size()
              << " queries=" << queries_
              << " lz_and_z_time=" << lz_and_z_time
              << " s=" << s_
              << " tau=" << tau_
              << " max_leaf_size=" << max_leaf_size_
              << " cut_first_levels=" << (cut_first_levels_ ? "true" : "false")
              << " extended_pruning=" << (extended_pruning_ ? "true" : "false")
              << " dynamic_prog=" << (dynamic_prog ? "true" : "false");


    auto const t1 = std::chrono::high_resolution_clock::now();
    auto* lpfPruned_bt_dp =
      new pasta::BlockTreeLPF<uint8_t, int32_t>(input_, tau_, max_leaf_size_, s_,
                                                    true /*always mark!*/,
                                                    cut_first_levels_,
						dynamic_prog, threads);

    auto const space_usage = lpfPruned_bt_dp->print_space_usage();
    
    auto const t2 = std::chrono::high_resolution_clock::now();


    lpfPruned_bt_dp->add_rank_support_omp(threads_);
    auto const space_usage_rs = lpfPruned_bt_dp->print_space_usage();
    
    auto const t3 = std::chrono::high_resolution_clock::now();


    auto const construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    auto const rank_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    std::cout << " space_usage=" << space_usage;
    std::cout << " construction_time_ms=" << construction_time;

    std::cout << " space_usage_with_rs=" <<  space_usage_rs;
    std::cout << " construction_time_with_rs_ms=" << construction_time + rank_construction_time;

    size_t result = 0;
    auto const t4 = std::chrono::high_resolution_clock::now();
    for (auto const& query : access_queries_) {
        result += lpfPruned_bt_dp->access(query);
    }
    auto const t5 = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < select_c_.size(); i++) {
        result += lpfPruned_bt_dp->rank(select_c_[i], access_queries_[i]);
    }

    auto const t6 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < select_c_.size(); i++) {
        result += lpfPruned_bt_dp->select(select_c_[i], select_queries_[i]);
    }

    auto const t7 = std::chrono::high_resolution_clock::now();


    auto const access_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count();
    auto const rank_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count();
    auto const select_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t6).count();
    
    std::cout << " access_time_ns=" << access_time
              << " rank_time_ns=" << rank_time
              << " select_time_ns=" << select_time
              << " threads=" << threads_
              << " check_sum=" << result << '\n';
}

  void run_construction_original() {

    std::string input_string(input_.begin(), input_.end());

    std::unordered_set<int> characters;
    for (char const c: input_) {
        characters.insert(c);
    }
    
    std::cout << "RESULT"
              << " algorithm=" << "belazzougui"
              << " input=" << input_path_
              << " prefix_size=" << input_.size()
              << " queries=" << queries_
              << " s=" << s_
              << " tau=" << tau_
              << " max_leaf_size=" << max_leaf_size_;


    auto const t1 = std::chrono::high_resolution_clock::now();

    auto* bt = new PBlockTree(input_string, tau_, max_leaf_size_);
    bt->process_back_pointers();
    bt->clean_unnecessary_expansions();
    bt->check();
    int64_t size = 0;
    for (auto const& level: bt->levelwise_iterator()) {
      size += static_cast<int64_t>(level.size()) * 12;
    }
    auto* abt = new PCBlockTree(bt);

    auto const t2 = std::chrono::high_resolution_clock::now();

    for (int c: characters) {
      bt->add_rank_select_support(c);
    }

    auto const t3 = std::chrono::high_resolution_clock::now();

    auto* cbt = new PCBlockTree(bt);
    
    auto const construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    auto const rank_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    std::cout << " space_usage=" << cbt_size_no_rank(*cbt);
    std::cout << " construction_time_ms=" << construction_time;

    std::cout << " space_usage_with_rs=" <<  size;
    std::cout << " construction_time_with_rs_ms=" << construction_time + rank_construction_time << '\n';
}

  
  void load_input() {
    std::ifstream stream(input_path_.c_str(), std::ios::in | std::ios::binary);
    if (!stream) {
      std::cerr << "File " << input_path_ << " not found\n";
      exit(1);
    }
    stream.seekg(0, std::ios::end);
    uint64_t size = stream.tellg();
    if (prefix_size_ > 0) {
      size = std::min(prefix_size_, size);
    }
    stream.seekg(0);
    input_.resize(size);
    stream.read(reinterpret_cast<char *>(input_.data()), size);
    stream.close();
  }

  void compute_number_lz_phrases() {
    if (number_lz_phrases_computed_) {
      return;
   }

    int32_t number_lz_phrases = 0;
    std::vector<int32_t> lpf2(input_.size());
    std::vector<int32_t> lpf_ptr2(input_.size());
    std::vector<int32_t> lz2;

    //lpf_array_omp(input_, lpf2, lpf_ptr2, threads_);
    pasta::lpf_array_ansv(input_,lpf2, lpf_ptr2, threads_);
    pasta::calculate_lz_factor(number_lz_phrases, lpf2, lz2);

    number_lz_phrases_ = number_lz_phrases;
    number_lz_phrases_computed_ = true;
  }

  void prepare_queries() {
    std::unordered_map<char, int64_t> hist;
    std::unordered_set<int> characters;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine(rnd_device());
    std::uniform_int_distribution<uint64_t> dist(0, input_.size() - 1);
    for (int i =  0; i < input_.size(); i++) {
        ++hist[input_[i]];
    }
    for (size_t i = 0; i < queries_; ++i) {
        access_queries_.push_back(dist(mersenne_engine));
    }
    for (size_t i = 0; i < queries_; ++i) {
        uint8_t x = 0;
        int xsum = 0;
        while (xsum + hist[x] < access_queries_[i]) {
            xsum += hist[x];
            x++;
        }
        if (access_queries_[i] - xsum > 0) {
            select_c_.push_back(x);
            select_queries_.push_back(access_queries_[i] - xsum);
        }
    }
  }

private:
  int cbt_size_no_rank(PCBlockTree& cbt) {
    int bt_bv_size = sizeof(void*);
    for (sdsl::bit_vector* bv : cbt.bt_bv_) {
        bt_bv_size += (int)sdsl::size_in_bytes(*bv);
    }
    int bt_bv_rank_size = sizeof(void*);
    for (sdsl::rank_support_v<1>* bvr : cbt.bt_bv_rank_) {
        bt_bv_rank_size +=(int)sdsl::size_in_bytes(*bvr);
    }
    int bt_offsets_size = sizeof(void*);
    for (sdsl::int_vector<>* offsets: cbt.bt_offsets_) {
        bt_offsets_size += (int)sdsl::size_in_bytes(*offsets);
    }
    int leaf_string_size = (int)sdsl::size_in_bytes(*cbt.leaf_string_);
    int alphabet_size = (int)sdsl::size_in_bytes(*cbt.alphabet_);
    int mapping_size = sizeof(int) * 256;
    // int partial_total_size = mapping_size + alphabet_size + bt_bv_size+ bt_bv_rank_size+ bt_offsets_size + leaf_string_size;
    return leaf_string_size + bt_offsets_size + bt_bv_rank_size + bt_bv_size + alphabet_size + mapping_size;
}

int cbt_size_rank(PCBlockTree& cbt) {
    int result = 0;
    int bt_prefix_ranks_first_level_size = (int)cbt.bt_first_level_prefix_ranks_.size()+1 * sizeof(void*);;
    for (auto pair: cbt.bt_first_level_prefix_ranks_) {
        bt_prefix_ranks_first_level_size += (int)sdsl::size_in_bytes(*(pair.second));
    }
    int bt_first_ranks_total_size = (int)(cbt.bt_first_ranks_.size()+1) * sizeof(void*);
    for (const auto& pair: cbt.bt_first_ranks_) {
        int size = 0;
        for (sdsl::int_vector<>* ranks: pair.second) {
            size += (int)sdsl::size_in_bytes(*ranks);
        }
        bt_first_ranks_total_size += size;
    }
    int bt_second_ranks_total_size = (cbt.bt_second_ranks_.size()+1) * sizeof(void*);
    for (auto pair: cbt.bt_second_ranks_) {
        int size = 0;
        for (sdsl::int_vector<>* ranks: pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_second_ranks_total_size += size;
    }
    int bt_prefix_ranks_total_size = (cbt.bt_prefix_ranks_.size()+1) * sizeof(void*);
    for (auto pair: cbt.bt_prefix_ranks_) {
        int size = 0;
        for (sdsl::int_vector<>* ranks: pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_prefix_ranks_first_level_size += size;
    }
    result = bt_prefix_ranks_first_level_size + bt_first_ranks_total_size + bt_second_ranks_total_size + bt_prefix_ranks_total_size;
    return result + cbt_size_no_rank(cbt);
}
  
}; // class BlockTreeBenchmark

int main(int argc, char* argv[]) {  
  tlx::CmdlineParser cp;

  cp.set_description("Block Tree Construction Benchmark");
  cp.set_author("Florian Kurpicz <florian@kurpicz.org>\n        Daniel Meyer");

  std::string input_path = "";
  cp.add_param_string("input", input_path, "Path to input file.");
  
  uint64_t size = 0;
  cp.add_bytes('p', "prefix", size, "Number of bytes of the prefix to process.");

  uint64_t queries = 1'000'000;
  cp.add_bytes('q', "queries", queries, "Number of queries tested.");
  
  uint64_t tau = 2;
  cp.add_bytes('t', "tau", tau, "Tuning parameter tau of the block tree, i.e., "
               "number of children of non-root and non-leaf nodes.");

  uint64_t s = 1;
  cp.add_bytes('s', "s", s, "Tuning parameter s of the block tree, i.e., number"
               " of children of the root node.");

  uint64_t max_leaf_size = 4;
  cp.add_bytes('m', "max_leaf_size", max_leaf_size, "Maximum size of a leaf "
               "before text is represented plainly.");
  
  bool cut_first_level = true;
  cp.add_flag('c', "cut_first_levels", cut_first_level , "Cut frist levels of "
              "block tree.");
    
  bool extended_pruning = true;
  cp.add_flag('e', "extended_pruning", extended_pruning, "Enable extended "
              "pruning.");

  uint64_t threads = 1;
  cp.add_bytes('h', "threads", threads, "Number of threads used.");

  bool original = false;
  cp.add_flag('o', "original", original, "Run original block tree "
	      "implementation only.");
  
  bool parallel = false;
  cp.add_flag('l', "parallel", parallel, "Run parallel implementation.");

  // process command line
  if (!cp.process(argc, argv)) {
    return -1; // some error occurred and help was always written to user.
  }

  BlockTreeBenchmark bench(input_path, size, queries, tau, s, max_leaf_size,
                           cut_first_level, extended_pruning, threads, original,
			   parallel);

  bench.run();

  return 0;
}
