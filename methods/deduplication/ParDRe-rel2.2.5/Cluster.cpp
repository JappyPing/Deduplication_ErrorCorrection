/*
 * Cluster.cpp
 *
 *  Created on: 16/06/2015
 *      Author: gonzales
 */

#include "Cluster.h"

Cluster::Cluster(uint8_t keyLength, uint8_t missMatch, int opticalDist){
	_keyLength = keyLength;
	_missMatch = missMatch;
	_onlyOptical = (opticalDist >= 0);
	_maxOpticalDist = opticalDist;
}

Cluster::~Cluster() {
	_seqs.clear();
}

void Cluster::includeSeq(Sequence *seq){

	_seqs.push_back(seq);
}

void Cluster::removeDup(){

	size_t numSeqs = _seqs.size();
	uint8_t *distTable = new uint8_t[numSeqs];
	uint8_t iterMis;

	Sequence *firstSeq = _seqs[0];

	for(size_t i=1; i<numSeqs; i++){
		iterMis = _dist(firstSeq, _seqs[i]);

		distTable[i] = iterMis;
		if(iterMis <= _missMatch){
			if(_equal(firstSeq, _seqs[i])){
				if(firstSeq->getAvgQual() >= _seqs[i]->getAvgQual()){
					_seqs[i]->_duplicate = true;
				} else {
					firstSeq->_duplicate = true;
				}
			}
		}
	}

	Sequence *it1, *it2;
	uint8_t diff;

	for(size_t i=1; i<numSeqs-1; i++){
		it1 = _seqs[i];
		if(!it1->_duplicate){
			for(size_t j=i+1; j<numSeqs; j++){
				it2 = _seqs[j];
				if(!it2->_duplicate){
					if(_onlyOptical && (it2->getOpticalCoord(true) - it1->getOpticalCoord(true) > _maxOpticalDist)){
						break;
					}

					diff = abs(distTable[i]-distTable[j]);
					if(diff<=_missMatch){
						if(_equal(it1, it2)){
							if(it1->getAvgQual() >= it2->getAvgQual()){
								it2->_duplicate = true;
							} else {
								it1->_duplicate = true;
							}
						}
					}
				}
			}
		}
	}

	delete [] distTable;
}


bool Cluster::_equal(Sequence *seq1, Sequence *seq2){

	if(seq1->_duplicate || seq2->_duplicate){
		return false;
	}

	if(seq1->_length != seq2->_length){
		return false;
	}

	if(_onlyOptical){
		if(!seq1->isEqualNameBase(seq2)){
			return false;
		}

		if(_calcOpticalDist(seq1, seq2) > _maxOpticalDist){
			return false;
		}
	}

	uint8_t numMissMatch = 0;
	uint64_t iter1, iter2, xorRes, popc;

	for(uint32_t i=0; i<seq1->_compressLength; i++){
		iter1 = seq1->_compressSuffix[i];
		iter2 = seq2->_compressSuffix[i];

		xorRes = iter1 ^ iter2;

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

uint8_t Cluster::_dist(Sequence *seq1, Sequence *seq2){

	if(seq1->_length != seq2->_length){
		return seq1->_length;
	}

	uint8_t numMissMatch = 0;
	uint64_t iter1, iter2, xorRes, popc;

	for(uint32_t i=0; i<seq1->_compressLength; i++){
		iter1 = seq1->_compressSuffix[i];
		iter2 = seq2->_compressSuffix[i];

		xorRes = iter1^iter2;

		popc = Utils::popcount(xorRes);

		numMissMatch += popc/2;
		if(popc%2){
			numMissMatch++;
		}
	}

	return numMissMatch;
}

int Cluster::_calcOpticalDist(Sequence *seq1, Sequence *seq2){
	int auxX = seq1->getOpticalCoord(true)-seq2->getOpticalCoord(true);
	auxX *= auxX;
	int auxY = seq1->getOpticalCoord(false)-seq2->getOpticalCoord(false);
	auxY *= auxY;

	return sqrt(auxX+auxY);
}
