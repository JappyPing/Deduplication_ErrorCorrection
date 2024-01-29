/*
 * Options.cpp
 *
 *  Created on: Jun 1, 2015
 *      Author: gonzalez
 */

#include "Options.h"
#include "Utils.h"

Options::Options() {
	/*set default values*/
	_setDefaults();
}
Options::~Options() {
	_readFileNames.clear();
	_outFileNames.clear();
}
void Options::_setDefaults() {

	_compressedFiles = 0;
	_readFileNames.clear();
	_outFileNames.clear();
	_numThreads = 1;
	_paired = false;

	_missMatch = 0;
	_numBasesCompared = -1;
	_prefixLen = 20;
	_print = true;
	_blockSize = 50000;
	_opticalDist = -1;
}

void Options::printUsage()
{
	fprintf(stderr, "\n%s version %s\n", PROGRAM_NAME, PROGRAM_VERSION);
	fprintf(stderr, "\nUsage: %s -i readFile1 [-p readFile2] [output] [options]\n", PROGRAM_NAME);
	/*the file input options*/
	fprintf(stderr, "Input:\n");
	fprintf(stderr,
			"\t-i <string> readFile1 (sequence file in FASTA/FASTQ format)\n");
	fprintf(stderr,
			"\t-p <string> readFile2 (sequence file in FASTA/FASTQ format for paired scenarios)\n");

	fprintf(stderr, "Output:\n");
	fprintf(stderr,
			"\t-o <string> outFile1 (output sequence file in FASTA/FASTQ format, default = readFile1.NonDup)\n");
	fprintf(stderr,
			"\t-r <string> outFile2 (second output sequence file in FASTA/FASTQ format for paired scenarios, default = readFile2.NonDup)\n");

	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-m <int> (number of allowed mismatches, default = %u)\n",
			_missMatch);
	fprintf(stderr, "\t-l <int> (prefix length used for computation, default = %u)\n",
			_prefixLen);
	fprintf(stderr, "\t-c <int> (number of bases to compare for each sequence (starting from the beginning), default is equal to the sequence length (all bases are compared))\n");
	fprintf(stderr, "\t-d <int> (in case you only want to discard optical duplicates, specify with this parameter the Euclidean distance that you consider two sequences are optical duplicates. By default, all duplicates or near-duplicates, either optical or PCR, are removed)\n");
	fprintf(stderr, "\t-b <int> (number of sequences that will be read in a block, default = %u)\n",
			_blockSize);
	fprintf(stderr, "\t-t <int> (number of threads per process, default = %d)\n",
			_numThreads);
	fprintf(stderr, "\t-z [<int>] (for gzip-compressed input and output, with optional compression level (1=fast, 9=best), if no level specified, uses gzip default; if 0 orno -z option, no compression.)\n");
	fprintf(stderr, "Others:\n");
	fprintf(stderr, "\t-np (do not write output files, just report on stderr)\n");
	fprintf(stderr, "\t-h  (print out the usage of the program)\n");
}

string Options::_makeOutputName(const string& inName)
{
    string aux = inName;
    if (_compressedFiles) {
        if (aux.substr(aux.length() - 3) == ".gz") {
            aux.resize(aux.length() - 3);
        }
    }
    aux.append(".NonDup_");
    aux.append(to_string(MPI::COMM_WORLD.Get_size()));
    if (_compressedFiles) {
        aux.append(".gz");
    }
    return (aux);
}

bool Options::parse(int argc, char* argv[])
{
	int intVal;
	int argind = 1;

	if (argc < 2) {
		Utils::log("Not enough parameters\n");
		return false;
	}

	/*check the help*/
	if (!strcmp(argv[argind], "-h") || !strcmp(argv[argind], "-?")) {
		return false;
	}

	/*print out the command line*/
	fprintf(stderr, "Command: ");
	for (int i = 0; i < argc; ++i) {
		fprintf(stderr, "%s ", argv[i]);
	}
	fputc('\n', stderr);

	bool withFile = false;
	bool withOutputFile = false;
	bool withPairedOutputFile = false;

	while (argind < argc) {
		/*single-end sequence files*/
		if (!strcmp(argv[argind], "-i")) {
			argind++;
			if (argind < argc) {
				_readFileNames.push_back(argv[argind]);
				argind++;
				withFile = true;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-p")) {
			argind++;
			if (argind < argc) {
				_readFileNames.push_back(argv[argind]);
				_paired = true;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-o")) {
			argind++;
			if (argind < argc) {
				string aux = argv[argind];
				_outFileNames.push_back(aux);
				withOutputFile = true;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-r")) {
			argind++;
			if (argind < argc) {
				string aux = argv[argind];
				_outFileNames.push_back(aux);
				withPairedOutputFile = true;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-m")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_missMatch = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-l")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				if (intVal > 27){
					Utils::log("WARNING: The length of the prefix can not be higher than 27. Using 27 instead");
					intVal = 27;
				}

				argind++;
				_prefixLen = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}

		} else if (!strcmp(argv[argind], "-c")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_numBasesCompared = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		}else if (!strcmp(argv[argind], "-d")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_opticalDist = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-b")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_blockSize = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-t")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_numThreads = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-z")) {
			intVal = 0;
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) == 1
						&& *(argv[argind]) >= '0'
						&& *(argv[argind]) <= '9') {
					sscanf(argv[argind], "%d", &intVal);
					argind++;
					_compressedFiles = intVal;
				}
                else {
                    _compressedFiles = Z_DEFAULT_COMPRESSION;
                }
			} else {
                _compressedFiles = Z_DEFAULT_COMPRESSION;
			}
		}
		else if (!strcmp(argv[argind], "-h")) {
			return false;
		} else if (!strcmp(argv[argind], "-np")) {
			_print = false;
			argind++;
		} else {
			Utils::log("Unknown parameter: %s\n", argv[argind]);
			return false;
		}
	}

	if(!withFile){
		Utils::log("No input file specified with -i\n");
		return false;
	}

	if(!withOutputFile){
		string aux = _makeOutputName(_readFileNames[0]);
		_outFileNames.push_back(aux);
	}

	if(!withPairedOutputFile && _paired){
		string aux = _makeOutputName(_readFileNames[1]);
		_outFileNames.push_back(aux);
	}

	if((_numBasesCompared >= 0) && (_numBasesCompared < _prefixLen)){
		_numBasesCompared = _prefixLen;
		Utils::log("The number of compared bases must be at least as high as the prefix length. New value: %d\n", _numBasesCompared);
	}

	return true;
}
