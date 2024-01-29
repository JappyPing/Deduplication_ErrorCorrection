/*
 * Sequence.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include "Utils.h"

struct Sequence
{
	/*member functions*/
	Sequence();
	Sequence(const Sequence& s);
	~Sequence();
	void setNameSize(size_t size);
	void setQualsSize(size_t size);
	inline void releaseQuals() {
		if (_quals) {
			delete[] _quals;
			_quals = NULL;
		}
	}

	inline uint64_t getKey(){
		return _key;
	}

	// This must be called when the name is set
	void setOpticalCoords();

	inline int getOpticalCoord(bool isX){
		if(isX){
			return _opticalX;
		} else {
			return _opticalY;
		}
	}

	inline bool isEqualNameBase(Sequence *seq2){
		for(uint32_t i=0; i<_nameSize; i++){
			if(_name[i] != seq2->_name[i]){
				return false;
			}

			if(_name[i] == '\0'){
				return true;
			}
		}

		return true;
	}

	uint64_t setKey(uint8_t *bases, uint8_t prefixLength, uint32_t length);
	void setSuffix(uint8_t *bases, int numBasesCompare);
	void decodeBases(uint8_t *bases);

	float getAvgQual();

	void clear();
	void printCompress(gzFile file, bool print_quals);
	void printUncompress(FILE * file, bool print_quals);

	/*member variables*/
	uint8_t* _name;
	uint8_t* _suffixName;
	uint8_t* _bases;
	uint8_t _prefixLength;
	uint64_t _key;
	uint64_t *_compressSuffix; // Only the number of bases that will be compared
	uint32_t _compressLength;
	uint8_t* _quals;
	uint32_t _length;
	uint32_t _nameSize;
	uint32_t _nameLength; // Really used
	uint32_t _suffixNameSize;
	int _opticalX;
	int _opticalY;
	bool _duplicate;

	float _avgQual;
};

struct SeqGreater{
	bool operator()(Sequence* lx, Sequence* rx ) {
		return (lx->getOpticalCoord(true) < rx->getOpticalCoord(true));
	}
};

#endif /* SEQUENCE_H_ */
