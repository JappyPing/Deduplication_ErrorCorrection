/*
 * PairedEnd.h
 *
 *  Created on: 22/06/2015
 *      Author: gonzales
 */

#ifndef PAIREDEND_H_
#define PAIREDEND_H_

#include "ClusterPaired.h"
#include "PairedSeqFileParser.h"

class PairedEnd {
public:
	PairedEnd(Options* options);
	virtual ~PairedEnd();

	virtual void execute();

private:
	uint32_t _getSeqs(Sequence* seqs1, Sequence* seqs2, int rank, int numP, uint8_t prefixLength);

protected:
	void _printOut(uint64_t myNumReads, double t);

	Options *_options;
	PairedSeqFileParser* _parser;
	SeqFileParser* _parserOut1;
	SeqFileParser* _parserOut2;
	uint64_t _numClusters;
	unordered_map<uint64_t, ClusterPaired*> _clHashTable;

	uint8_t _prefixLength;
	int _numBasesComp;
	uint8_t _missMatch;
	uint32_t _blockSize;
	int _opticalDist;

	// To access the map
	int _numTh;
	mutex _mutex;
	bool _withLock;
};

#endif /* PAIREDEND_H_ */
