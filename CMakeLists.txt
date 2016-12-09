cmake_minimum_required(VERSION 3.6)
project(SSJOINLib)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER  gcc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lboost_program_options -lboost_system -std=c++11 -pthread -static -Wwrite-strings")


set(SOURCE_FILES ssjoin-0.1/cmdline_mpjoin_def.h ssjoin-0.1/cmdline_groupjoin_dice.cxx ssjoin-0.1/cmdline_mpjoin_hamming.cxx ssjoin-0.1/cmdline_mpjoin_jaccard.cxx ssjoin-0.1/adaptjoin.h ssjoin-0.1/adaptjoin_impl.h ssjoin-0.1/adaptjoin_indexes.h ssjoin-0.1/allpairs.h ssjoin-0.1/allpairs_policies.h ssjoin-0.1/candidateset.h ssjoin-0.1/classes.h ssjoin-0.1/cmdline.h ssjoin-0.1/cmdline_adaptjoin.cxx ssjoin-0.1/cmdline_allpairs.cxx ssjoin-0.1/cmdline_groupjoin.cxx ssjoin-0.1/cmdline_groupjoin_cosine.cxx ssjoin-0.1/cmdline_groupjoin_def.h ssjoin-0.1/cmdline_groupjoin_dice.cxx ssjoin-0.1/cmdline_groupjoin_hamming.cxx ssjoin-0.1/cmdline_groupjoin_jaccard.cxx ssjoin-0.1/cmdline_mpjoin.h ssjoin-0.1/cmdline_mpjoin.cxx ssjoin-0.1/cmdline_groupjoin_jaccard.cxx ssjoin-0.1/cmdline_mpjoin_cosine.cxx ssjoin-0.1/cmdline_mpjoin_dice.cxx
        ssjoin-0.1/config.h.in ssjoin-0.1/cpucycles.h ssjoin-0.1/data.h ssjoin-0.1/frequencysorting.h ssjoin-0.1/functions.h ssjoin-0.1/GitSHA1.cxx.in ssjoin-0.1/GitSHA1.h ssjoin-0.1/groupjoin.h ssjoin-0.1/groupjoin_indexes.h ssjoin-0.1/indexes.h ssjoin-0.1/input.cxx ssjoin-0.1/input.h ssjoin-0.1/intify.cxx
        ssjoin-0.1/inv_index.h ssjoin-0.1/lengthfilter.h ssjoin-0.1/lsh.cxx ssjoin-0.1/lsh.h ssjoin-0.1/mpjoin.h ssjoin-0.1/mpjoinpolicies.h ssjoin-0.1/mpltricks.h ssjoin-0.1/output.h ssjoin-0.1/ppjoinpolicies.h ssjoin-0.1/similarity.h ssjoin-0.1/statistics.cxx ssjoin-0.1/statistics.h ssjoin-0.1/template_unroll.h ssjoin-0.1/timing.cxx ssjoin-0.1/timing.h ssjoin-0.1/verify.h main.cpp)
add_executable(SSJOINLib ${SOURCE_FILES})