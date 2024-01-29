/*
 * ClusterPired.cpp
 *
 *  Created on: 22/06/2015
 *      Author: gonzales
 */

#include "ClusterPaired.h"

ClusterPaired::ClusterPaired(uint8_t keyLength, uint8_t missMatch, int opticalDist){
	_keyLength = keyLength;
	_missMatch = missMatch;
	_onlyOptical = (opticalDist >= 0);
	_maxOpticalDist = opticalDist;
}

ClusterPaired::~ClusterPaired() {

	_seqs1.clear();
	_seqs2.clear();
}

void ClusterPaired::includeSeq(Sequence* seq1, Sequence* seq2){

	_seqs1.push_back(seq1);
	_seqs2.push_back(seq2);
}

void ClusterPaired::removeDup(){

	size_t numSeqs = _seqs1.size();
	uint8_t *distTable = new uint8_t[numSeqs];
	uint8_t iterMis;

	Sequence *firstSeq1 = _seqs1[0];
	Sequence *firstSeq2 = _seqs2[0];

	for(size_t i=1; i<numSeqs; i++){
		iterMis = _dist(firstSeq1, firstSeq2, _seqs1[i], _seqs2[i]);
		distTable[i] = iterMis;

		if(iterMis <= _missMatch){
			if(_equal(firstSeq1, firstSeq2, _seqs1[i], _seqs2[i])){
				if((firstSeq1->getAvgQual()+firstSeq2->getAvgQual()) >= (_seqs1[i]->getAvgQual()+_seqs2[i]->getAvgQual())){
					_seqs1[i]->_duplicate = true;
					_seqs2[i]->_duplicate = true;
				} else {
					firstSeq1->_duplicate = true;
					firstSeq2->_duplicate = true;
				}
			}
		}
	}

	Sequence *itOrigin1, *itOrigin2, *itNew1, *itNew2;
	uint8_t diff;

	for(size_t i=1; i<numSeqs-1; i++){
		itOrigin1 = _seqs1[i];
		itOrigin2 = _seqs2[i];
		if(!itOrigin1->_duplicate){
			for(size_t j=i+1; j<numSeqs; j++){
				itNew1 = _seqs1[j];
				itNew2 = _seqs2[j];
				if(!itNew1->_duplicate){
					if(_onlyOptical && (itNew1->getOpticalCoord(true) - itOrigin1->getOpticalCoord(true) > _maxOpticalDist)){
						break;
					}

					diff = abs(distTable[i]-distTable[j]);
					if(diff<=_missMatch){
						if(_equal(itOrigin1, itOrigin2, itNew1, itNew2)){
							if((itOrigin1->getAvgQual()+itOrigin2->getAvgQual()) >= itNew1->getAvgQual()+itNew2->getAvgQual()){
								itNew1->_duplicate = true;
								itNew2->_duplicate = true;
							} else {
								itOrigin1->_duplicate = true;
								itOrigin2->_duplicate = true;
							}
						}
					}
				}
			}
		}
	}

	delete [] distTable;
}

bool ClusterPaired::_equal(Sequence* seqOrigin1, Sequence* seqOrigin2, Sequence* seqNew1, Sequence* seqNew2){

	if(seqOrigin1->_duplicate || seqOrigin2->_duplicate || seqNew1->_duplicate || seqNew2->_duplicate){
		return false;
	}

	if((seqOrigin1->_length != seqNew1->_length) || (seqOrigin2->_length != seqNew2->_length)){
		return false;
	}

	if(_onlyOptical){
		if(!seqOrigin1->isEqualNameBase(seqNew1)){
			return false;
		}

		if(_calcOpticalDist(seqOrigin1, seqNew1) > _maxOpticalDist){
			return false;
		}
	}

	uint8_t numMissMatch = 0;
	uint64_t iterOrigin, iterNew, xorRes, popc;

	for(uint32_t i=0; i<seqOrigin1->_compressLength; i++){
		iterOrigin = seqOrigin1->_compressSuffix[i];
		iterNew = seqNew1->_compressSuffix[i];

		xorRes = iterOrigin ^ iterNew;
		popc = Utils::popcount(xorRes);
		numMissMatch += popc/2;
		if(popc%2){
			numMissMatch++;
		}

		if(numMissMatch > _missMatch){
			return false;
		}
	}

	// The paired sequence
	for(uint32_t i=0; i<seqOrigin2->_compressLength; i++){
		iterOrigin = seqOrigin2->_compressSuffix[i];
		iterNew = seqNew2->_compressSuffix[i];

		xorRes = iterOrigin ^ iterNew;
		popc = Utils::popcount(xorRes);
		numMissMatch += popc/2;
		if(popc%2){
			numMissMatch++;
		}

		if(numMissMatch > _missMatch){
			return false;
		}
	}

	return true;
}

uint8_t ClusterPaired::_dist(Sequence* seqOrigin1, Sequence* seqOrigin2, Sequence* seqNew1, Sequence* seqNew2){

	if(seqOrigin1->_duplicate || seqOrigin2->_duplicate || seqNew1->_duplicate || seqNew2->_duplicate){
		return false;
	}

	if((seqOrigin1->_length != seqNew1->_length) || (seqOrigin2->_length != seqNew2->_length)){
		return false;
	}

	uint8_t numMissMatch = 0;
	uint64_t iterOrigin, iterNew, xorRes, popc;

	for(uint32_t i=0; i<seqOrigin1->_compressLength; i++){
		iterOrigin = seqOrigin1->_compressSuffix[i];
		iterNew = seqNew1->_compressSuffix[i];

		xorRes = iterOrigin ^ iterNew;
		popc = Utils::popcount(xorRes);
		numMissMatch += popc/2;
		if(popc%2){
			numMissMatch++;
		}
	}

	// The paired sequence
	for(uint32_t i=0; i<seqOrigin2->_compressLength; i++){
		iterOrigin = seqOrigin2->_compressSuffix[i];
		iterNew = seqNew2->_compressSuffix[i];

		xorRes = iterOrigin ^ iterNew;
		popc = Utils::popcount(xorRes);
		numMissMatch += popc/2;
		if(popc%2){
			numMissMatch++;
		}
	}

	return numMissMatch;
}

int ClusterPaired::_calcOpticalDist(Sequence *seq1, Sequence *seq2){
	int auxX = seq1->getOpticalCoord(true)-seq2->getOpticalCoord(true);
	auxX *= auxX;
	int auxY = seq1->getOpticalCoord(false)-seq2->getOpticalCoord(false);
	auxY *= auxY;

	return sqrt(auxX+auxY);
}
