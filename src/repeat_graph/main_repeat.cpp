//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>

#include "../sequence/sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include "../repeat_graph/repeat_graph.h"
#include "../repeat_graph/multiplicity_inferer.h"
#include "../repeat_graph/haplotype_resolver.h"
#include "../repeat_graph/graph_processing.h"
#include "../repeat_graph/repeat_resolver.h"
#include "../repeat_graph/output_generator.h"
#include "../contigger/contig_extender.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inAssembly, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads, 
			   std::string& configPath, bool& unevenCov,
			   std::string& matchMode, float& ovlpDivergence, bool& keepHaplotypes)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: flye-repeat "
				  << " --input-seq path --out-dir path --config path\n"
				  << "\t\t[--reads pathi] [--log path] [--treads num] [--kmer size] \n"
				  << "\t\t[--min-ovlp size] [--max-divergence X] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --input-asm path\tpath to input sequences\n"
				  << "  --out-dir path\tpath to output dir\n"
				  << "  --config path\tpath to the config file\n\n"
				  << "Optional arguments:\n"
				  << "  --reads reads to simplify the graph \n"
				  << "  --kmer size\tk-mer size [default = 15] \n"
				  << "  --min-ovlp size\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "  --debug \t\tenable debug output "
				  << "[default = false] \n"
				  << "  --meta \t\tenable uneven coverage (metagenome) mode "
				  << "[default = false] \n"
				  << "  --keep-haplotypes \t\tdo not collapse alternative haplotypes "
				  << "[default = false] \n"
				  << "  --log log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "  --match-mode\teither of [local, semi, dovetail] "
				  << "[default = local] \n"
				  << "  --ovlp-divergence max overlap divergence "
				  << "[default = 0.05] \n"
				  << "  --threads num_threads\tnumber of parallel threads \n";
	};
	
	int optionIndex = 0;
	static option longOptions[] =
	{
		{"input-seq", required_argument, 0, 0},
		{"reads", required_argument, 0, 0},
		{"out-dir", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"max-divergence", required_argument, 0, 0},
		{"match-mode", required_argument, 0, 0},
		{"meta", no_argument, 0, 0},
		{"keep-haplotypes", no_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{0, 0, 0, 0}
	};

	int opt = 0;
	while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex)) != -1)
	{
		switch(opt)
		{
		case 0:
			if (!strcmp(longOptions[optionIndex].name, "kmer"))
				kmerSize = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "threads"))
				numThreads = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-ovlp"))
				minOverlap = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "match-mode"))
				matchMode = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "max-divergence"))
				ovlpDivergence = atof(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "log"))
				logFile = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "debug"))
				debug = true;
			else if (!strcmp(longOptions[optionIndex].name, "input-seq"))
				inAssembly = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "meta"))
				unevenCov = true;
			else if (!strcmp(longOptions[optionIndex].name, "keep-haplotypes"))
				keepHaplotypes = true;
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-dir"))
				outFolder = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "config"))
				configPath = optarg;
			break;

		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (outFolder.empty() || inAssembly.empty() || configPath.empty())
	{
		printUsage();
		return false;
	}

	return true;
}

int repeat_main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	bool debugging = false;
	size_t numThreads = 1;
	int kmerSize = 15;
	int minOverlap = 1000;
	bool isMeta = true;
	bool keepHaplotypes = false; 
	std::string readsFasta;
	std::string inAssembly;
	std::string outFolder;
	std::string logFile;
	std::string configPath;
	std::string matchModeStr = "local";
	float maxOvlpDivergence = 0.05f;
	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inAssembly,
				   kmerSize, minOverlap, debugging, 
				   numThreads, configPath, isMeta, 
				   matchModeStr, maxOvlpDivergence, keepHaplotypes))  return 1;
	
	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);
	
	Logger::get().debug() << "Total RAM: " 
		<< getMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Available RAM: " 
		<< getFreeMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Total CPUs: " << std::thread::hardware_concurrency();

	Config::load(configPath);
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Parameters::get().minimumOverlap = minOverlap;
	Parameters::get().unevenCoverage = isMeta;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().info() << "Selected minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[isMeta];

	Logger::get().info() << "Reading sequences";

	SequenceContainer seqAssembly; 
	SequenceContainer seqReads;
	std::vector<std::string> readsList = splitString(readsFasta, ',');
	try
	{
		seqAssembly.loadFromFile(inAssembly);
		for (auto& readsFile : readsList)
		{
			seqReads.loadFromFile(readsFile);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	seqAssembly.buildPositionIndex();
	seqReads.buildPositionIndex();

	SequenceContainer edgeSequences;
	RepeatGraph rg(seqAssembly, &edgeSequences);
	GraphProcessor proc(rg, seqAssembly);
	ReadAligner aligner(rg, seqReads);
	OutputGenerator outGen(rg);

	Logger::get().info() << "Building repeat graph";

	OverlapDetector::MatchMode matchMode(OverlapDetector::MatchLocal);
	if (matchModeStr == "semi") matchMode = OverlapDetector::MatchSemiDovetail;
	if (matchModeStr == "dovetail") matchMode = OverlapDetector::MatchDovetail;
	Logger::get().info() << "Matching mode: " << matchModeStr;
	Logger::get().info() << "Max overlap divergence: " << maxOvlpDivergence;
	rg.build(matchMode, maxOvlpDivergence);
	rg.updateEdgeSequences();
	//outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_raw.gv");
	
	if (!readsList.empty())
	{
		Logger::get().info() << "Mapping reads on the graph";
		aligner.alignReads();
		MultiplicityInferer multInf(rg, aligner, seqAssembly);
		multInf.estimateCoverage();
		//multInf.splitNodes();

		RepeatResolver repResolver(rg, seqAssembly, seqReads, aligner, multInf);
		repResolver.resolveSimpleRepeats();
		//repResolver.findRepeats();
	}

	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/repeat_graph.gv");
	outGen.outputGfa(proc.getEdgesPaths(), outFolder + "/repeat_graph.gfa");
	//outGen.outputFasta(proc.getEdgesPaths(), outFolder + "/graph_edges.fasta");

	SequenceContainer emptyContainer;
	ContigExtender extender(rg, aligner, emptyContainer, emptyContainer);
	extender.generateUnbranchingPaths();

	outGen.outputDot(extender.getUnbranchingPaths(),
					 outFolder + "/repeat_graph_compact.gv");
	//outGen.outputFasta(extender.getUnbranchingPaths(),
	//				   outFolder + "/graph_edge_compact.fasta");
	outGen.outputGfa(extender.getUnbranchingPaths(),
					 outFolder + "/repeat_graph_compact.gfa");

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
