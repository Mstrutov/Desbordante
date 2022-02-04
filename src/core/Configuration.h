#pragma once

#include <string>

struct Configuration {
    bool is_find_keys = true;
    bool is_find_fds = true;
    bool is_null_equal_null = true;

    std::string ucc_error_measure = "g1prime";

    //Error settings
    double error_dev = 0;
    bool is_estimate_only = false;
    double max_ucc_error = 0.01;          // both for FD and UCC actually

    //Traversal settings
    int parallelism = 0;
    int max_threads_per_search_space = -1;
    bool is_defer_failed_launch_pads = true;
    std::string launch_pad_order = "error";

    unsigned int max_lhs = -1;


    //Sampling settings
    unsigned int sample_size = 10000;
    double sample_booster = 10;
    int seed = 0; //check seed, mb use boost::optional
    double estimate_confidence = 0; //       -||-
    int random_ascend_threads = 2;

    //Cache settings
    double caching_probability = 0.5;
    unsigned int nary_intersection_size = 4;

    //Miscellaneous settings
    bool is_check_estimates = false;
    bool is_initial_pause = false;
    std::string fd_score_measure = "none";
};