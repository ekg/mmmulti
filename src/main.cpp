#include <iostream>
#include <random>
#include <limits>
#include <cassert>
#include <chrono>
#include "args.hxx"
#include <utility>
#include "mmmultimap.hpp"
#include "mmmultiset.hpp"
#include "mmiitree.hpp"
#include "paryfor.hpp"

int main(int argc, char** argv) {

    args::ArgumentParser parser("memmapped multimap interface");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    //args::ValueFlag<std::string> in_file(parser, "FILE", "use this input file for a uint64_t sort", {'i', "in"});
    args::ValueFlag<std::string> test_file(parser, "FILE", "test mmmultimap with random data in this file", {'T', "test-file"});
    args::ValueFlag<uint64_t> test_size(parser, "N", "test this many pairs", {'s', "test-size"});
    args::ValueFlag<uint64_t> max_val(parser, "N", "generate test data in the range [1,max_value]", {'M', "max-value"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    args::ValueFlag<uint64_t> unique_value_tests(parser, "N", "number of unique value calls to make", {'u', "unique-vals"});
    args::Flag test_multiset(parser, "multiset", "test the multiset", {'m', "test-multiset"});
    args::Flag test_complex(parser, "complex", "test the multimap with complex values", {'c', "test-complex-values"});
    args::Flag test_unpadded(parser, "unpadded", "test the multimap without padding for random access", {'P', "test-unpadded"});
    args::Flag test_iitree(parser, "iitree", "test the iitree from cranges", {'I', "test-iitree"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    int num_threads = 1;
    if (args::get(threads)) {
        num_threads = args::get(threads);
    }

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // these both stack smash! why??
    //std::mt19937 gen(100); //Standard mersenne_twister_engine seeded with rd()
    //std::mt19937_64 gen(100); //Standard mersenne_twister_engine seeded with rd();
    std::default_random_engine gen(rd());

    if (!args::get(test_file).empty() && !args::get(test_multiset) && !args::get(test_iitree)) {
        if (args::get(test_complex)) {
            uint64_t max_key = args::get(max_val);
            std::uniform_int_distribution<uint64_t> dis(1, max_key);
            //std::vector<uint64_t> x; x.reserve(1e8);
            std::remove(args::get(test_file).c_str());
            auto mm = std::make_unique<mmmulti::map<uint64_t, std::pair<uint64_t, uint64_t>>>(args::get(test_file), std::make_pair(0,0));
            uint64_t x_len = args::get(test_size);
            mm->open_writer();
            paryfor::parallel_for<uint64_t>(
                0, x_len, num_threads, 1000,
                [&](uint64_t n) {
                    mm->append(dis(gen), std::make_pair(dis(gen), dis(gen)));
                });
            mm->index(num_threads, max_key);
            uint64_t i = 0;
            uint64_t key_count = 0;
            uint64_t value_count = 0;
            uint64_t unique_value_count = 0;
            bool first = true;
            uint64_t last = 0;
            mm->for_each_pair([&](const uint64_t& a, const std::pair<uint64_t, uint64_t>& b) {
                    if (first || a > last) {
                        ++key_count;
                        last = a;
                        first = false;
                        mm->for_unique_values_of(a, [&](const std::pair<uint64_t, uint64_t>& v) {
                                ++unique_value_count;
                            });
                    }
                    ++value_count;
                });
            // exercise unique value search
            uint64_t unique_value_test_count = args::get(unique_value_tests);
            auto start = std::chrono::system_clock::now();
            uint64_t x = 0;
            for (uint64_t i = 0; i < unique_value_test_count; ++i) {
                uint64_t q = dis(gen);
                mm->for_unique_values_of(q, [&](const std::pair<uint64_t, uint64_t>& v) {
                        ++x;
                    });
            }
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;

            std::cerr << key_count << " keys" << std::endl;
            std::cerr << value_count << " values" << std::endl;
            std::cerr << unique_value_count << " unique pairs" << std::endl;
            if (unique_value_test_count) std::cerr << elapsed_seconds.count()/unique_value_test_count << "s per unique value call" << std::endl;
            mm.reset();

        }
        else if (args::get(test_unpadded)) {
            uint64_t max_key = args::get(max_val);
            std::uniform_int_distribution<uint64_t> dis(1, max_key);
            //std::vector<uint64_t> x; x.reserve(1e8);
            std::remove(args::get(test_file).c_str());
            auto mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(args::get(test_file), 0);
            uint64_t x_len = args::get(test_size);
            mm->open_writer();
            paryfor::parallel_for<uint64_t>(
                0, x_len, num_threads, 1000,
                [&](uint64_t n) {
                    mm->append(dis(gen), dis(gen));
                });
            mm->index(num_threads);
            uint64_t i = 0;
            uint64_t key_count = 0;
            uint64_t value_count = 0;
            uint64_t unique_value_count = 0;
            bool first = true;
            uint64_t last = 0;
            mm->for_each_pair([&](const uint64_t& a, const uint64_t& b) {
                    if (first || a > last) {
                        ++key_count;
                        last = a;
                        first = false;
                        mm->for_unique_values_of(a, [&](const uint64_t& v) {
                                ++unique_value_count;
                            });
                    }
                    ++value_count;
                });
            std::cerr << key_count << " keys" << std::endl;
            std::cerr << value_count << " values" << std::endl;
            std::cerr << unique_value_count << " unique pairs (these can't be found with unpadded mmmultimaps)" << std::endl;
        }
        else {
            uint64_t max_key = args::get(max_val);
            std::uniform_int_distribution<uint64_t> dis(1, max_key);
            //std::vector<uint64_t> x; x.reserve(1e8);
            std::remove(args::get(test_file).c_str());
            auto mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(args::get(test_file), 0);
            uint64_t x_len = args::get(test_size);
            mm->open_writer();
            paryfor::parallel_for<uint64_t>(
                0, x_len, num_threads, 1000,
                [&](uint64_t n) {
                    mm->append(dis(gen), dis(gen));
                });
            mm->index(num_threads, max_key);
            uint64_t i = 0;
            uint64_t key_count = 0;
            uint64_t value_count = 0;
            uint64_t unique_value_count = 0;
            bool first = true;
            uint64_t last = 0;
            mm->for_each_pair([&](const uint64_t& a, const uint64_t& b) {
                    if (first || a > last) {
                        ++key_count;
                        last = a;
                        first = false;
                        mm->for_unique_values_of(a, [&](const uint64_t& v) {
                                ++unique_value_count;
                            });
                    }
                    ++value_count;
                });
            // exercise unique value search
            uint64_t unique_value_test_count = args::get(unique_value_tests);
            auto start = std::chrono::system_clock::now();
            uint64_t x = 0;
            for (uint64_t i = 0; i < unique_value_test_count; ++i) {
                uint64_t q = dis(gen);
                mm->for_unique_values_of(q, [&](const uint64_t& v) {
                        ++x;
                    });
            }
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;

            std::cerr << key_count << " keys" << std::endl;
            std::cerr << value_count << " values" << std::endl;
            std::cerr << unique_value_count << " unique pairs" << std::endl;
            if (unique_value_test_count) std::cerr << elapsed_seconds.count()/unique_value_test_count << "s per unique value call" << std::endl;
        }
    } else if (!args::get(test_file).empty() && args::get(test_multiset)) {
        uint64_t max_value = args::get(max_val);
        std::uniform_int_distribution<uint64_t> dis(1, max_value);
        //std::vector<uint64_t> x; x.reserve(1e8);
        std::remove(args::get(test_file).c_str());
        auto ms = std::make_unique<mmmulti::set<uint64_t>>(args::get(test_file));
        uint64_t x_len = args::get(test_size);
        ms->open_writer();
        paryfor::parallel_for<uint64_t>(
            0, x_len, num_threads, 1000,
            [&](uint64_t n) {
                //Use dis to transform the random unsigned int generated by gen into an int in the range
                ms->append(dis(gen));
            });
        ms->index(num_threads);
        uint64_t i = 0;
        uint64_t value_count = 0;
        uint64_t unique_value_count = 0;
        uint64_t sum1 = 0;
        // exercise unique value search
        ms->for_each_value_count([&](const uint64_t& value, const uint64_t& count) {
                ++unique_value_count;
                value_count += count;
                sum1 += count * value;
            });
        uint64_t second_count = 0;
        uint64_t sum2 = 0;
        for (auto v = ms->begin(); v != ms->end(); ++v) {
            ++second_count;
            sum2 += *v;
        }
        std::cerr << value_count << " values, expected " << x_len << std::endl;
        std::cerr << unique_value_count << " unique pairs" << std::endl;
        std::cerr << "sums " << sum1 << " " << sum2 << std::endl;
    } else if (!args::get(test_file).empty() && args::get(test_iitree)) {
        //} else if (args::get(test_iitree)) {
        std::remove(args::get(test_file).c_str());
        auto tree = std::make_unique<mmmulti::iitree<uint64_t, uint64_t>>(args::get(test_file));
        uint64_t max_value = args::get(max_val);
        std::uniform_int_distribution<uint64_t> dis(0, max_value);
        uint64_t x_len = args::get(test_size);
        tree->open_writer();
        paryfor::parallel_for<uint64_t>(
            0, x_len, num_threads, 1000,
            [&](uint64_t n) {
                uint64_t q = dis(gen);
                uint64_t r = std::min(q + 100, max_value);
                tree->add(q, r, dis(gen));
            });
        tree->index(num_threads);
        for (int n=0; n<max_value; ++n) {
            if (n % 1000 == 0) std::cerr << n << "\r";
            tree->overlap(
                n, n+1,
                [&](const uint64_t& start,
                    const uint64_t& end,
                    const uint64_t& data) {
                    if (start > n || end < n) {
                        std::cerr << "tree broken at " << n << std::endl;
                    }
                });
        }
        std::cerr << max_value << std::endl;
        //tree->for_each_entry([&](const mmmulti::iitree<uint64_t, uint64_t>::Interval& ival) {
        //        std::cerr << ival.st << ".." << ival.en << " " << ival.data << std::endl;
        //    });
        //std::cerr << std::endl;
    }

    return 0;

}
