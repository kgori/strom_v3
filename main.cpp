#include <iostream>
#include "node.hpp"
#include "tree.hpp"
#include "libhmsbeagle/beagle.h"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    std::cout << "Starting..." << std::endl;
    Tree tree;
    std::cout << "\nFinished!" << std::endl;

    int ntaxa = 4;
    bool use_tip_partials = true;
    int nsites = 100;
    int mtrxCount = 1;
    int nrates = 1;
    int *rsrcList = new int[1];
    int rsrcCnt = 1;
    bool sse_vectorization = true;
    bool auto_scaling = true;
    long requirementFlags = 0;
    BeagleInstanceDetails instDetails;
    beagleCreateInstance(
            ntaxa,		// tipCount
            ntaxa + 2,	// partialsBufferCount
            (use_tip_partials ? 0 : ntaxa),			// compactBufferCount
            4, 			// stateCount
            nsites,		// patternCount
            1,			// eigenBufferCount
            mtrxCount,	// matrixBufferCount,
            nrates,     // categoryCount
            3,          // scalingBuffersCount
            rsrcList,	// resourceList
            rsrcCnt,	// resourceCount
            (sse_vectorization ? BEAGLE_FLAG_VECTOR_SSE : 0) | (auto_scaling ? BEAGLE_FLAG_SCALING_AUTO : 0),         // preferenceFlags
            requirementFlags,			// requirementFlags
            &instDetails);

    return 0;
}