#pragma once
#include <string>

struct Configuration {
    //Error settings
    double errorDev = 0;
    bool isEstimateOnly = false;

    //Traversal settings
    int parallelism = 0;
    int maxThreadsPerSearchSpace = -1;
    bool isDeferFailedLaunchPads = true;
    std::string launchPadOrder = "error";

    //Sampling settings
    int sampleSize = 10000;
    double sampleBooster = 10;
    unsigned int seed = 0; //check seed, mb use boost::optional
    double estimateConfidence = 0; //       -||-
    int randomAscendThreads = 2;

    //Cache settings
    double cachingProbability = 0.5;
    int naryIntersectionSize = 4;

    //Miscellaneous settings
    bool isCheckEstimates = false;
    bool isInitialPause = false;
    std::string fdScoreMeasure = "none";
};