/*
 * Cluster.h
 *
 *  Created on: 16/06/2015
 *      Author: gonzales
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "Sequence.h"
#include "Options.h"

class Cluster {
public:
	Cluster(uint8_t keyLength, uint8_t missMatch, int opticalDist);
	virtual ~Cluster();

	void includeSeq(Sequence *seq);

	void removeDup();

	inline int numSeqs(){
		return _seqs.size();
	}

	inline Sequence* getSeq(int pos){
		return _seqs[pos];
	}

	inline void sort(){
		std::sort(_seqs.begin(), _seqs.end(), SeqGreater());
	}

private:
	// Indicate if two sequences are equal for clustering
	bool _equal(Sequence *seq1, Sequence *seq2);

	// Indicate the distance between two sequences
	uint8_t _dist(Sequence *seq1, Sequence *seq2);

	int _calcOpticalDist(Sequence *seq1, Sequence *seq2);

	uint8_t _keyLength;
	uint8_t _missMatch;
	vector<Sequence*> _seqs;

	bool _onlyOptical;
	int _maxOpticalDist;
};

#endif /* CLUSTER_H_ */
