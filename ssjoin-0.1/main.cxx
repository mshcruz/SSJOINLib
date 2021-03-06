/* Copyright 2014-2015 Willi Mann
 *
 * This file is part of set_sim_join.
 *
 * set_sim_join is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with set_sim_join.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include "functions.h"
#include "classes.h"
#include "cmdline.h"
#include "GitSHA1.h"

#include "input.h"

//algorithms
//#include "ppjoin.h"
//#include "ppjoinpolicies.h"
//#include "bitjoin.h"
#include "lsh.h"


#include "timing.h"

#include "output.h"


//typedef PPJoinAndPlus<Jaccard, IndexOnTheFlyPolicy, PPJoinPolicy> PPJoin;
//typedef PPJoinAndPlus<Jaccard, IndexOnTheFlyPolicy, PPJoinPlusPolicy> PPJoinPlus;

typedef LSH<Jaccard> JLSH;



namespace po = boost::program_options;

int main(int argc, char ** argv) {
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("timings", "output timings")
		("statistics", "output statistics")
		("mpjoin", "Use mpjoin filtering")
		("pljoin", "Use pljoin filtering - length filtering based on lookup position")
		("indexfirst", "index first, then join")
		("result", po::value<std::string>(), "print result to file specified")
		("algorithm", po::value<std::string>(), "algorithm to test")
		("threshold", po::value<double>(), "jaccard threshold")
		("input", po::value<std::string>(), "file, each line a record")
		("foreign-input", po::value<std::string>(), "file, each line a record for foreign join")
		("whitespace", "preprocess input: split by whitespace")
		("qgram", po::value<unsigned int>(), "preprocess input: create q-grams of length q")
		("cosine", "use cosine similarity")
		("dice", "use dice similarity")
		("hamming", "use hamming similarity")
		("jaccard", "use jaccard similarity (default)")
		("allprefext", "Support all prefix extensions (AdaptJoin)")
		("suffixfilter", "suffix filter pre-verification")
		("version", "output version information (git SHA1)")
		;

	po::positional_options_description p;
	p.add("input", 1);
	p.add("algorithm", 1);
	p.add("threshold", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
          options(desc).positional(p).run(), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 0;
	}

	if(vm.count("version")) {
		std::cout << g_GIT_SHA1;
#ifdef NO_STAT_COUNTERS
		std::cout << "-NO_STAT_COUNTERS";
#endif
#ifdef CAND_ONLY
		std::cout << "-CAND_ONLY";
#endif
#ifdef CYCLE_COUNT
		std::cout << "-CYCLE_COUNT";
#endif
#ifdef CYCLE_COUNT_SF
		std::cout << "-CYCLE_COUNT_SF";
#endif
#ifdef EXT_STATISTICS
		std::cout << "-EXT_STATISTICS";
#endif
#ifdef LONG_VERIFICATION
		std::cout << "-LONG_VERIFICATION";
#endif
#ifdef MPJ_LINKEDLIST
		std::cout << "-MPJ_LINKEDLIST";
#endif
		std::cout << std::endl;
		return 0;
	}

	if(!vm.count("threshold")) {
		std::cerr << "no threshold specified";
		std::cout << desc << "\n";
	}

	double threshold = vm["threshold"].as<double>();
	
	if (!vm.count("algorithm")) {
		std::cout << "algorithm was not specified " << std::endl;
		std::cout << desc << "\n";
		return 1;
	}
	std::string algorithm = vm["algorithm"].as<std::string>();

	Algorithm * algo = NULL;

	if (algorithm == "ppjoin") {
		algo = mpjoin_cmd_line(vm);
	} else if (algorithm == "groupjoin") {
		algo = groupjoin_cmd_line(vm);
	} else if (algorithm == "allpairs") {
		algo = allpairs_cmd_line(vm);
	} else if (algorithm == "adaptjoin") {
		algo = adaptjoin_cmd_line(vm);
#ifdef LSHAVAILABLE
	} else if (algorithm == "lsh") {
		if(vm.count("cosine")) {
			std::cout << "lsh does not support cosine similarity" << std::endl;
			exit(1);
		} else if (vm.count("dice")) {
			std::cout << "lsh does not support dice similarity" << std::endl;
			exit(1);
		}
		algo = new JLSH(threshold);
#endif
	} else {
		std::cout << "Algorithm " << algorithm << " not available or unknown" << std::endl;
		exit(1);
	}
	
	if (!vm.count("threshold") ) {
		std::cout << "threshold was not specified " << std::endl;
		std::cout << desc << "\n";
		return 1;
	} 
	if (!vm.count("input")) {
		std::cout << "inputfile was not specified " << std::endl;
		std::cout << desc << "\n";
		return 1;
	} 

	unsigned int qgrams = 0;
	bool rawnumbers = true;
	if(vm.count("qgram") || vm.count("whitespace")) {
		rawnumbers = false;

		if(vm.count("qgram")) {
			qgrams = vm["qgram"].as<unsigned int>();
		}
	}

	IntRecords records;
	IntRecords foreignrecords;
	bool foreign = false;
	
	std::string inputfile = vm["input"].as<std::string>();

	Timing timings;

	if(rawnumbers) {
		std::cout << "Reading Input ..." << std::endl;
		Timing::Interval * readrawrecords = timings.add("readrawrecords");
		add_raw_input<algoaddrecord>().add_input(algo, inputfile, algoaddrecord());
		timings.finish(readrawrecords);

		if(vm.count("foreign-input")) {
			Timing::Interval * readrawforeignrecords = timings.add("readrawforeignrecords");

			std::string foreigninputfile = vm["foreign-input"].as<std::string>();
			add_raw_input<algoaddforeignrecord>().add_input(algo, foreigninputfile, algoaddforeignrecord());

			timings.finish(readrawforeignrecords);
		}
	} else {
		//First, common preparations

		std::cout << "Preprocessing Input ..." << std::endl;

		Timing::Interval * commonprepare = timings.add("commonprepare");
		get_int_records_fast gir(algo);
		// Fill the records DS with integer records that are prepared for the prefix filtering
		if(qgrams == 0) {
			tokenize_whitespace tw;
			gir.start<tokenize_whitespace>(inputfile, tw);
		} else {
			tokenize_qgram tq(qgrams);
			gir.start<tokenize_qgram>(inputfile, tq);
		}
		if(vm.count("foreign-input")) {
			foreign = true;
			std::string foreignfilename = vm["foreign-input"].as<std::string>();
			if(qgrams == 0) {
				tokenize_whitespace tw;
				gir.foreignsets<tokenize_whitespace>(foreignfilename, tw);
			} else {
				tokenize_qgram tq(qgrams);
				gir.foreignsets<tokenize_qgram>(foreignfilename, tq);
			}

		}
		timings.finish(commonprepare);

		// Second, algorithm specific preparations

		Timing::Interval * algopreparerecords = timings.add("algoPrepRecs");
		algo->preparerecords();
		timings.finish(algopreparerecords);

		if(foreign) {
			Timing::Interval * algoprepareforeignrecords = timings.add("algoPrepForRecs");
			algo->prepareforeignrecords();
			timings.finish(algoprepareforeignrecords);
		}

		Timing::Interval * algopreparefinished = timings.add("algoPrepCleanUp");
		algo->preparefinished();
		timings.finish(algopreparefinished);
	}

	// Prepare output handler
	HandleOutput * handleoutput;
	HandleOutputPairs * handleoutputpairs = NULL;

	if (vm.count("result")) {
		handleoutputpairs = new HandleOutputPairs();
		handleoutput = handleoutputpairs;
	} else {
		handleoutput = new HandleOutputCount();
	}

	std::cout << "running algorithm " << algorithm << " with threshold " << threshold << std::endl;

	Timing::Interval * algoindex = timings.add("algoindex");
	algo->doindex();
	timings.finish(algoindex);

	Timing::Interval * algojoin = timings.add("algojoin");
	algo->dojoin(handleoutput);
	timings.finish(algojoin);


	if(vm.count("result")) {
		assert(handleoutputpairs != NULL);
		std::string resultfile = vm["result"].as<std::string>();
		print_result<print_pairs_lf, HandleOutputPairs::output_type>(resultfile, handleoutputpairs->getOutput(), print_pairs_lf());
	}

	std::cout << "Result pair count: " << handleoutput->getCount() << std::endl;

	if(vm.count("statistics")) {
		std::cout << algo->statistics;
		std::cout << extStatistics;
	}

	if(algo != NULL) {
		delete algo;
	}

	delete handleoutput;

	if(vm.count("timings")) {
		std::cout << "Timings:" << std::endl;
		std::cout << timings;
	}
	return 0;
}
