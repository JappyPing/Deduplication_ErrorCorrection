/*
 * ClusterPaired.h
 *
 *  Created on: 22/06/2015
 *      Author: gonzales
 */

#ifndef CLUSTERPAIRED_H_
#define CLUSTERPAIRED_H_

#include "Sequence.h"
#include "Options.h"

class ClusterPaired {
public:
	ClusterPaired(uint8_t keyLength, uint8_t missMatch, int opticalDist);
	virtual ~ClusterPaired();

	void includeSeq(Sequence* seq1, Sequence* seq2);

	void removeDup();

	inline int numSeqs(){
		return _seqs1.size();
	}

	inline void getSeq(int pos, Sequence *&seq1, Sequence *&seq2){

		seq1 = _seqs1[pos];
		seq2 = _seqs2[pos];
	}

	inline void sort(){
		std::sort(_seqs1.begin(), _seqs1.end(), SeqGreater());
		std::sort(_seqs2.begin(), _seqs2.end(), SeqGreater());
	}

private:
	// Indicate if two pairs of sequences are equal for clustering
	bool _equal(Sequence* seqOrigin1, Sequence* seqOrigin2, Sequence* seqNew1, Sequence* seqNew2);

	// Indicate the distance between two sequences
	uint8_t _dist(Sequence* seqOrigin1, Sequence* seqOrigin2, Sequence* seqNew1, Sequence* seqNew2);

	int _calcOpticalDist(Sequence *seq1, Sequence *seq2);

	uint8_t _keyLength;
	uint8_t _missMatch;
	vector<Sequence*> _seqs1;
	vector<Sequence*> _seqs2;

	bool _onlyOptical;
	int _maxOpticalDist;
};

#endif /* CLUSTERPAIRED_H_ */
