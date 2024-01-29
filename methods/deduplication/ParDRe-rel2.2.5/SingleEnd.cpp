/*
 * SingleEnd.cpp
 *
 *  Created on: 16/06/2015
 *      Author: gonzales
 */

#include "SingleEnd.h"

SingleEnd::SingleEnd(Options* options) {

	int myRank = MPI::COMM_WORLD.Get_rank();

	_options = options;
	_prefixLength = options->getPrefixLen();
	_numBasesComp = options->getNumBasesCompared();
	_missMatch = options->getMissMatch();
	_blockSize = options->getBlockSize();
	_opticalDist = options->getOpticalDist();

	_numTh = options->getNumThreads();
	_withLock = _numTh>1;

	_parser = new SeqFileParser(options, options->getInputFileList()[0].c_str(), false);

	if(_options->getPrint()){
		string outPath = options->getOutFileList()[0];
		if(!myRank){
			Utils::log("Process %d/%d: To print in %s\n", myRank, MPI::COMM_WORLD.Get_size(), outPath.c_str());
		}
	/*if(myRank){
		outPath.append("_id");
		outPath.append(to_string(myRank));
	}
	_parserOut = new SeqFileParser(options, outPath.c_str(), true, false);
	_parserOut->setFormat(_parser->getFormat());*/
	}
}

SingleEnd::~SingleEnd() {
	delete _parser;
	//delete _parserOut;
	_clHashTable.clear();
}

uint32_t SingleEnd::_getSeqs(Sequence* seqs, int rank, int numP, uint8_t prefixLength){

	uint32_t numReads = 0;

	while(numReads < _blockSize){
		int seqLength = _parser->getSeq(seqs[numReads], numP, rank, prefixLength, _numBasesComp);

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

void _funcRemoveDup(unordered_map<uint64_t, Cluster*>::iterator* it, unordered_map<uint64_t, Cluster*>::iterator end,
		mutex *m, bool withLock, int id){

	Cluster* cl;

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

void SingleEnd::execute(){

	double stime;
	stime = Utils::getSysTime();

    int myRank = MPI::COMM_WORLD.Get_rank();
    int numP = MPI::COMM_WORLD.Get_size();

	Sequence* seqs = new Sequence[_blockSize];
    Sequence* iterSeq;
	uint32_t readsBlock;
	uint64_t myNumReads = 0;
	uint64_t key;
	unordered_map<uint64_t, Cluster*>::iterator it;
	Cluster* cl;

	while((readsBlock = _getSeqs(seqs, myRank, numP, _prefixLength)) > 0){
		for(uint32_t i=0; i<readsBlock; i++){
			iterSeq = &seqs[i];
			key = iterSeq->getKey();

			it = _clHashTable.find(key);
			if(it == _clHashTable.end()){
				cl = new Cluster(_prefixLength, _missMatch, _opticalDist);
				_clHashTable.emplace(key, cl);
			} else {
				cl = it->second;
			}

			cl->includeSeq(iterSeq);
		}

		seqs = new Sequence[_blockSize];
		myNumReads += readsBlock;

		if(!(myNumReads % 1000000)){
			Utils::log("Process %d/%d: %lu clustered reads in %.2f seconds\n",
					myRank, numP, myNumReads, Utils::getSysTime()-stime);
		}
	}

	if(_opticalDist > 0){
		for(auto it = _clHashTable.begin(); it != _clHashTable.end(); it++){
			cl = it->second;
			cl->sort();
		}
	}

	Utils::log("Process %d/%d: finished the clustering of %lu reads in %.2f seconds\n",
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

	Utils::log("Process %d/%d: Time to analyze %d reads with %d threads: %.2f seconds\n",
			myRank, numP, myNumReads, _numTh, t);

	_printOut(myNumReads, t);
}

void SingleEnd::_printOut(uint64_t myNumReads, double t){
    // Print the information of the clusters
    uint64_t printReads = 0;
    Cluster *cl;
    Sequence *iterSeq;

    int myRank = MPI::COMM_WORLD.Get_rank();
    int numP = MPI::COMM_WORLD.Get_size();

    // Not all processes can print at the same time
    for(int i=0; i<numP; i++){
    	if(myRank == i){
    		string outPath = _options->getOutFileList()[0];
    		_parserOut = new SeqFileParser(_options, outPath.c_str(), true, i, false);

    		_parserOut->setFormat(_parser->getFormat());

    		for(auto it = _clHashTable.begin(); it != _clHashTable.end(); it++){
    			cl = it->second;
    			for(int s=0; s<cl->numSeqs(); s++){
    				iterSeq = cl->getSeq(s);
    				if(!iterSeq->_duplicate){
    					if(_options->getPrint()){
    						_parserOut->printSeq(iterSeq);
    					}
    					printReads++;
    				}
    			}
    			delete cl;
    		}

			_parserOut->flushOutput();
			delete _parserOut;
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
