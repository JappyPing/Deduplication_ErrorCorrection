/*
 * SeqFileParser.h
 *
 *  Created on: Aug 13, 2015
 *      Author: gonzalez
 */

#ifndef PAIREDSEQFILEPARSER_H_
#define PAIREDSEQFILEPARSER_H_

#include "SeqFileParser.h"

class PairedSeqFileParser
{
public:
	/*public member functions*/
	PairedSeqFileParser(Options* options, const char* path1, const char* path2, size_t BUFFER_SIZE = 4095);
	~PairedSeqFileParser();

	inline int getFormat(){
		return _parser1->getFormat();
	}

	inline int getSeq(Sequence &seq1, Sequence &seq2, int numProcs, int myId, uint8_t prefixLength, int numBasesComp){

		int length1 =_parser1->getSeq(seq1, numProcs, myId, prefixLength, numBasesComp);

		/*int secondNumBases;
		if(numBasesComp < 0){
			secondNumBases = numBasesComp;
		} else if(numBasesComp-length1 < 0){
			secondNumBases = 0;
		} else {
			secondNumBases = numBasesComp-length1;
		}*/

		int length2 = _parser2->getSeq(seq2, 1, 0, 0, numBasesComp);

		if(length1 == -1){
			return -1;
		}

		if(!length1 || !length2){
			return 0;
		}

		return length1;
	}

private:
	SeqFileParser* _parser1;
	SeqFileParser* _parser2;
};

#endif /* PAIREDSEQFILEPARSER_H_ */
