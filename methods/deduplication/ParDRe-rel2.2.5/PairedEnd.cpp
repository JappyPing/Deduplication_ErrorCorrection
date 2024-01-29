/*
 * PairedEnd.cpp
 *
 *  Created on: 22/06/2015
 *      Author: gonzales
 */

#include "PairedEnd.h"

PairedEnd::PairedEnd(Options* options) {
	_options = options;
	_numTh = options->getNumThreads();
	_prefixLength = options->getPrefixLen();
	_numBasesComp = options->getNumBasesCompared();
	_missMatch = options->getMissMatch();
	_blockSize = options->getBlockSize();
	_opticalDist = options->getOpticalDist();

	_numTh = options->getNumThreads();
	_withLock = _numTh>1;

	_parser = new PairedSeqFileParser(options, options->getInputFileList()[0].c_str(),
			options->getInputFileList()[1].c_str());

	int myRank = MPI::COMM_WORLD.Get_rank();

	if(!myRank){
		Utils::log("To print in %s and %s\n", options->getOutFileList()[0].c_str(), options->getOutFileList()[1].c_str());
	}

	/*if(_options->getPrint()){
	string outPath = options->getOutFileList()[0];
	if(myRank){
		outPath.append("_id");
		outPath.append(to_string(myRank));
	}
	_parserOut1 = new SeqFileParser(options, outPath.c_str(), true, false);
	_parserOut1->setFormat(_parser->getFormat());

	outPath = options->getOutFileList()[1];
	if(myRank){
		outPath.append("_id");
		outPath.append(to_string(myRank));
	}
	_parserOut2 = new SeqFileParser(options, outPath.c_str(), true, false);
	_parserOut2->setFormat(_parser->getFormat());
	}*/
}

PairedEnd::~PairedEnd() {
	delete _parser;
	//delete _parserOut1;
	//delete _parserOut2;
	_clHashTable.clear();
}

uint32_t PairedEnd::_getSeqs(Sequence* seqs1, Sequence* seqs2, int rank, int numP, uint8_t prefixLength){

	uint32_t numReads = 0;

	while(numReads < _blockSize){
		int seqLength = _parser->getSeq(seqs1[numReads], seqs2[numReads], numP, rank, prefixLength, _numBasesComp);

		if(seqLength == 0){
			break;
		}

		if(seqLength == -1){
			continue;
		}

		numReads++;
	}

	return numReads;
}

void _funcRemoveDup(unordered_map<uint64_t, ClusterPaired*>::iterator* it,
		unordered_map<uint64_t, ClusterPaired*>::iterator end,
		mutex *m, bool withLock, int id){

	ClusterPaired* cl;

	while(1){
		if(withLock){
			m->lock();
		}

		if((*it) == end){
			if(withLock){
				m->unlock();
			}
			break;
		}

		cl = (*it)->second;
		(*it)++;

		if(withLock){
			m->unlock();
		}

		cl->removeDup();
	}
}

void PairedEnd::execute(){

	double stime = Utils::getSysTime();

    int myRank = MPI::COMM_WORLD.Get_rank();
    int numP = MPI::COMM_WORLD.Get_size();

	Sequence* seqs1 = new Sequence[_blockSize];
	Sequence* seqs2 = new Sequence[_blockSize];
    Sequence* iterSeq1;
    Sequence* iterSeq2;
	uint32_t readsBlock;
	uint64_t myNumReads = 0;
	uint64_t key;
	unordered_map<uint64_t, ClusterPaired*>::iterator it;
	ClusterPaired* cl;

	while((readsBlock = _getSeqs(seqs1, seqs2, myRank, numP, _prefixLength)) > 0){
		for(uint32_t i=0; i<readsBlock; i++){
			iterSeq1 = &seqs1[i];
			iterSeq2 = &seqs2[i];
			key = iterSeq1->getKey();

			it = _clHashTable.find(key);
			if(it == _clHashTable.end()){
				cl = new ClusterPaired(_prefixLength, _missMatch, _opticalDist);
				_clHashTable.emplace(key, cl);
			} else {
				cl = it->second;
			}

			cl->includeSeq(iterSeq1, iterSeq2);
		}

		seqs1 = new Sequence[_blockSize];
		seqs2 = new Sequence[_blockSize];
		myNumReads += readsBlock;

		if(!(myNumReads % 1000000)){
			Utils::log("Process %d/%d: %lu clustered paired reads in %.2f seconds\n",
					myRank, numP, myNumReads, Utils::getSysTime()-stime);
		}
	}

	if(_opticalDist > 0){
		for(auto it = _clHashTable.begin(); it != _clHashTable.end(); it++){
			cl = it->second;
			cl->sort();
		}
	}

	Utils::log("Process %d/%d: finished the clustering of %lu paired reads in %.2f seconds\n",
			myRank, numP, myNumReads, Utils::getSysTime()-stime);

    vector<thread> threads;
    it = _clHashTable.begin();

    for(int i=0; i<_numTh; i++){
        threads.push_back(thread(_funcRemoveDup, &it, _clHashTable.end(),
        		&_mutex, _withLock, i));
    }

    for(int i=0; i<_numTh; i++){
        threads[i].join();
    }

	double t = Utils::getSysTime()-stime;

	Utils::log("Process %d/%d: Time to analyze %d paired reads with %d threads: %.2f seconds\n",
			myRank, numP, myNumReads, _numTh, t);

	_printOut(myNumReads, t);
}

void PairedEnd::_printOut(uint64_t myNumReads, double t){

    // Print the information of the clusters
    ClusterPaired *cl;
    Sequence *seq1;
    Sequence *seq2;
    uint64_t printReads = 0;

    int myRank = MPI::COMM_WORLD.Get_rank();
    int numP = MPI::COMM_WORLD.Get_size();

    for(int i=0; i<numP; i++){
    	if(myRank == i){
    		string outPath = _options->getOutFileList()[0];
    		_parserOut1 = new SeqFileParser(_options, outPath.c_str(), true, i, false);
    		_parserOut1->setFormat(_parser->getFormat());

    		outPath = _options->getOutFileList()[1];
    		_parserOut2 = new SeqFileParser(_options, outPath.c_str(), true, i, false);
    		_parserOut2->setFormat(_parser->getFormat());

    		// Not all processes can print at the same time
    		for(auto it = _clHashTable.begin(); it != _clHashTable.end(); it++){
    			cl = it->second;
    			for(int s=0; s<cl->numSeqs(); s++){
    				cl->getSeq(s, seq1, seq2);
    				if(!seq1->_duplicate){
    					if(_options->getPrint()){
    						_parserOut1->printSeq(seq1);
    						_parserOut2->printSeq(seq2);
    					}
    					printReads++;
    				}
    			}
    			delete cl;
    		}
			_parserOut1->flushOutput();
			_parserOut2->flushOutput();
			//fflush(NULL);
    	}
    	MPI::COMM_WORLD.Barrier();
    }

   	/*if(_options->getPrint()){
   	MPI::COMM_WORLD.Barrier();
    stime = Utils::getSysTime();
   	if(!myRank){
   		for(int i=1; i<numP; i++){
			string outputFile = _options->getOutFileList()[0];
			outputFile.append("_id");
			outputFile.append(to_string(i));
			string command = "cat ";
			command.append(outputFile);
			command.append(" >> ");
			command.append(_options->getOutFileList()[0]);
			system(command.c_str());

			command = "rm ";
			command.append(outputFile);
			system(command.c_str());

			outputFile = _options->getOutFileList()[1];
			outputFile.append("_id");
			outputFile.append(to_string(i));
			command = "cat ";
			command.append(outputFile);
			command.append(" >> ");
			command.append(_options->getOutFileList()[1]);
			system(command.c_str());

			command = "rm ";
			command.append(outputFile);
			system(command.c_str());
    	}

   		Utils::log("SUMMARY: Print time %.2f\n", printTime+Utils::getSysTime()-stime);
    }}*/

   	// To gather information from all processes
   	uint64_t totalReads, totalPrints, minReads, maxReads;
   	double minTime, maxTime;

   	MPI::COMM_WORLD.Reduce(&myNumReads, &totalReads, 1, MPI::UNSIGNED_LONG, MPI::SUM, 0);
   	MPI::COMM_WORLD.Reduce(&myNumReads, &minReads, 1, MPI::UNSIGNED_LONG, MPI::MIN, 0);
   	MPI::COMM_WORLD.Reduce(&myNumReads, &maxReads, 1, MPI::UNSIGNED_LONG, MPI::MAX, 0);
   	MPI::COMM_WORLD.Reduce(&printReads, &totalPrints, 1, MPI::UNSIGNED_LONG, MPI::SUM, 0);
   	MPI::COMM_WORLD.Reduce(&t, &minTime, 1, MPI::DOUBLE, MPI::MIN, 0);
   	MPI::COMM_WORLD.Reduce(&t, &maxTime, 1, MPI::DOUBLE, MPI::MAX, 0);

   	if(!myRank){
   		float percent = (float) totalPrints;
   		percent /= (float) totalReads;
   		percent *= 100.00;

   		Utils::log("SUMMARY: Non duplicated paired reads %lu/%lu (%.2f %%)\n", totalPrints, totalReads, percent);

   		percent = (float) maxReads;
   		percent -= (float) minReads;
   		percent /= (float) maxReads;
   		percent *= 100.00;

   		Utils::log("SUMMARY: Min/Max analyzed reads per process: %lu/%lu (%.2f %% of imbalance)\n",
   				minReads, maxReads, percent);

   		percent = (float) maxTime;
   		percent -= (float) minTime;
   		percent /= (float) maxTime;
   		percent *= 100.00;

   		Utils::log("SUMMARY: Min/Max runtime per process: %.2f/%.2f (%.2f %% of imbalance)\n",
   				minTime, maxTime, percent);
   	}
}
