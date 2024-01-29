/*
 * Sequence.cpp
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#include "Sequence.h"
#include "Utils.h"

const uint8_t decode[9] = { 'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T' };

const uint8_t code[26] = { 1, 0, 2, 0, 0, 0, 4, //A -> G
		0, 0, 0, 0, 0, 0, 0, //H->N
		0, 0, 0, 0, 0, 8, 0, //O->U
		0, 0, 0, 0, 0 //V->Z
		};

const uint8_t decodePrefix[5] = { 'A', 'C', 'G', 'T', 'N' };

const uint8_t codePrefix[26] = { 0, 4, 1, 4, 4, 4, 2, //A -> G
		4, 4, 4, 4, 4, 4, 4, //H->N
		4, 4, 4, 4, 4, 3, 4, //O->U
		4, 4, 4, 4, 4 //V->Z
		};


Sequence::Sequence() {
	_name = NULL;
	_suffixName = NULL;
	_bases = NULL;
	_quals = NULL;
	_key = 0;
	_prefixLength = 0;
	_compressSuffix = NULL;
	_compressLength = 0;
	_nameSize = 0;
	_nameLength = 0;
	_suffixNameSize = 0;
	_length = 0;
	_avgQual = -1;
	_duplicate = false;
	_opticalX = -1;
	_opticalY = -1;
}
Sequence::Sequence(const Sequence & s) {
	_length = s._length;
	_nameSize = s._nameSize;
	_nameLength = s._nameLength;
	_suffixNameSize = s._suffixNameSize;
	_compressLength = s._compressLength;
	_key = s._key;
	_prefixLength = s._prefixLength;
	_duplicate = s._duplicate;
	_avgQual = s._avgQual;
	_opticalX = s._opticalX;
	_opticalY = s._opticalY;
	if (_length == 0) {
		_name = NULL;
		_compressSuffix = NULL;
		_quals = NULL;
		return;
	}
	if (s._name) {
		_name = new uint8_t[_nameSize];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		strcpy((char*) _name, (const char*) s._name);
	}
	if (s._suffixName) {
		_suffixName = new uint8_t[_suffixNameSize];
		if (_suffixName == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		strcpy((char*) _suffixName, (const char*) s._suffixName);
	}
	if (s._compressLength) {
		_compressSuffix = new uint64_t[_compressLength];
		if (_compressSuffix == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_compressSuffix, s._compressSuffix, _compressLength*sizeof(uint64_t));
	}
	if (s._bases) {
		_bases = new uint8_t[_length];
		if (_bases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_bases, s._bases, _length);
	}
	if (s._quals) {
		_quals = new uint8_t[_length];
		if (_quals == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_quals, s._quals, _length);
	}
}
Sequence::~Sequence() {
	clear();
}
void Sequence::clear() {
	if (_name) {
		delete[] _name;
	}
	if (_suffixName) {
		delete[] _suffixName;
	}
	if (_compressSuffix) {
		delete[] _compressSuffix;
	}
	if (_quals != NULL) {
		delete[] _quals;
	}
	if (_bases != NULL) {
		delete[] _bases;
	}

	_name = NULL;
	_suffixName = NULL;
	_quals = NULL;
	_compressSuffix = NULL;
	_length = 0;
	_prefixLength = 0;
	_nameSize = 0;
	_nameLength = 0;
	_suffixNameSize = 0;
	_key = 0;
	_avgQual = -1;
	_duplicate = false;
}

uint64_t Sequence::setKey(uint8_t *bases, uint8_t prefixLength, uint32_t length){
	_key = 0;
	_length = length;
	int limit = prefixLength;
	if(prefixLength > length){
		limit = length;
	}

	_prefixLength = limit;

	for(int i=0; i<limit; i++){
		_key *= 5;
		_key += codePrefix[bases[i]];
	}

	return _key;
}

void Sequence::setSuffix(uint8_t *bases, int numBasesCompare){

	_bases = new uint8_t[_length];
	if (_bases == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	memcpy(_bases, bases, _length);

	uint8_t *iterBases = &bases[_prefixLength];

	// 16 bases per element
	uint32_t suffixLength;
	if(numBasesCompare >= 0){
		suffixLength = numBasesCompare-_prefixLength;
	} else {
		suffixLength = _length-_prefixLength;
	}

	_compressLength = suffixLength/16;
	if(suffixLength%16){
		_compressLength++;
	}

	if(_compressSuffix){
		delete [] _compressSuffix;
	}

	_compressSuffix = (uint64_t *) malloc(_compressLength*sizeof(uint64_t));
	for(uint32_t i=0; i<suffixLength/16; i++){
		_compressSuffix[i] = 0;
		for(uint32_t j=0; j<16; j++){
			_compressSuffix[i] <<= 4;
			_compressSuffix[i] += code[*iterBases];
			iterBases++;
		}
	}

	// The last block can be no complete
	if(suffixLength%16){
		_compressSuffix[_compressLength-1] = 0;
		for(uint32_t i=0; i<suffixLength%16; i++){
			_compressSuffix[_compressLength-1] <<= 4;
			_compressSuffix[_compressLength-1] += code[*(iterBases+i)];
		}
	}

	return;

}

void Sequence::decodeBases(uint8_t *bases){
	uint8_t *iterBases = &bases[_prefixLength-1];
	uint64_t auxKey = _key;
	for(uint32_t i=0; i<_prefixLength; i++){
		*iterBases = decodePrefix[auxKey%5];
		auxKey /= 5;
		iterBases--;
	}

	// decode the suffix
	iterBases = &bases[_length-1];
	uint64_t *iterSuffix = &_compressSuffix[_compressLength-1];
	uint32_t suffixLength = _length-_prefixLength;
	uint32_t numDecoded = 0;

	if(suffixLength%16){
		for(uint32_t i=0; i<suffixLength%16; i++){
			*iterBases = decode[(*iterSuffix)&15];
			*iterSuffix >>= 4;
			iterBases--;
		}
		iterSuffix--;
		numDecoded++;
	}

	for(; numDecoded < _compressLength; numDecoded++){
		for(uint32_t j=0; j<16; j++){
			*iterBases = decode[(*iterSuffix)&15];
			*iterSuffix >>= 4;
			iterBases--;
		}
		iterSuffix--;
	}
}

float Sequence::getAvgQual(){

	if((_avgQual==-1) && (_quals != NULL)){
		_avgQual = 0;
		for(uint32_t i=0; i<_length; i++){
			_avgQual += _quals[i];
		}
		_avgQual /= (float) _length;
	}

	return _avgQual;
}

void Sequence::setNameSize(size_t size) {
	_nameLength = size;
	if (size >= _nameSize) {
		_nameSize = size * 2;
		if (_name) {
			delete[] _name;
		}
		_name = new uint8_t[_nameSize];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
}
void Sequence::setQualsSize(size_t size) {

	if(_quals){
		delete [] _quals;
	}

	_quals = new uint8_t[size];
	if (_quals == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
}

void Sequence::printCompress(gzFile file, bool print_quals) {
	//print the sequence name
	/*if (_quals) {
		fputc('@', file);
	} else {
		fputc('>', file);
	}*/

	gzprintf(file, "%s", _name);
	if(_opticalX > 0){
		gzprintf(file, ":%d:%d%s", _opticalX, _opticalY, _suffixName);
	}
	gzprintf(file, "\n");

	//print the query sequence
	for (uint32_t i = 0; i < _length; ++i) {
		gzputc(file, decode[code[_bases[i]]]);
	}
	gzputc(file, '\n');

	//print the quality scores if available
	if (print_quals) {
		gzputc(file, '+');
		gzputc(file, '\n');
		for (uint32_t i = 0; i < _length; ++i) {
			gzputc(file, _quals[i]);
		}
		gzputc(file, '\n');
	}
}

void Sequence::printUncompress(FILE *file, bool print_quals) {
	//print the sequence name
	/*if (_quals) {
		fputc('@', file);
	} else {
		fputc('>', file);
	}*/

	fprintf(file, "%s", _name);
	if(_opticalX > 0){
		fprintf(file, ":%d:%d%s", _opticalX, _opticalY, _suffixName);
	}
	fprintf(file, "\n");

	//print the query sequence
	for (uint32_t i = 0; i < _length; ++i) {
		fputc(decode[code[_bases[i]]], file);
	}
	fputc('\n', file);

	//print the quality scores if available
	if (print_quals) {
		fputc('+', file);
		fputc('\n', file);
		for (uint32_t i = 0; i < _length; ++i) {
			fputc(_quals[i], file);
		}
		fputc('\n', file);
	}
}


void Sequence::setOpticalCoords(){
	int lastDigit = 0;
	// first scan forward to find end of first "word" on defline
	while (!isspace(_name[lastDigit]) && _name[lastDigit])
		lastDigit++;
	// then back up to first colon
	int firstDigit = lastDigit-1;
	while (_name[firstDigit] != ':' && firstDigit>0){
		if(!isdigit(_name[firstDigit])){
			lastDigit = firstDigit;
		}
		firstDigit--;
	}
	if (firstDigit == 0) {
		Utils::exit("Name '%s' does not have :X:Y coordinates\n", _name);
	}
	firstDigit++;

	_suffixNameSize = _nameLength-lastDigit-1;

	_suffixName = new uint8_t[_suffixNameSize+1];
	if (_suffixName == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	memcpy(_suffixName, &_name[lastDigit], _suffixNameSize*sizeof(uint8_t));
	_suffixName[_suffixNameSize] = '\0';

	_opticalY = 0;

	for(int i=firstDigit; i<lastDigit; i++){
		_opticalY *= 10;
		_opticalY += _name[i]-48;
	}

	// unsure: this was firstDigit-2, but should be -1?
	lastDigit = firstDigit-1;
	firstDigit = lastDigit-1;

	while(_name[firstDigit] != ':' && firstDigit>0){
		firstDigit--;
	}
	if (firstDigit == 0) {
		Utils::exit("Name '%s' does not have :X:Y coordinates\n", _name);
	}
	firstDigit++;

	_opticalX = 0;

	for(int i=firstDigit; i<lastDigit; i++){
		_opticalX *= 10;
		_opticalX += _name[i]-48;
	}

	// trim the name to its prefix (before the X and Y coordinate)
	_name[firstDigit-1] = '\0';
}
