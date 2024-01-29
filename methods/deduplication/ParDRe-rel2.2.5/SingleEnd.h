/*
 * SingleEnd.h
 *
 *  Created on: 16/06/2015
 *      Author: gonzales
 */

#ifndef SINGLEEND_H_
#define SINGLEEND_H_

#include "Cluster.h"
#include "SeqFileParser.h"

#include <fstream>

class SingleEnd {
public:
	SingleEnd(Options* options);
	virtual ~SingleEnd();

	virtual void execute();

private:
	uint32_t _getSeqs(Sequence* seqs, int rank, int numP, uint8_t prefixLength);

protected:

	void _printOut(uint64_t myNumReads, double t);

	Options *_options;
	SeqFileParser* _parser;
	SeqFileParser* _parserOut;
	uint64_t _numClusters;
	unordered_map<uint64_t, Cluster*> _clHashTable;

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

#endif /* SINGLEEND_H_ */
