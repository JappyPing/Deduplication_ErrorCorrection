/*
 * SeqFileParser.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef SEQFILEPARSER_H_
#define SEQFILEPARSER_H_

#include "Macros.h"
#include "Sequence.h"
#include "Utils.h"
#include "Options.h"

class SeqFileParser
{
public:
	/*public member functions*/
	SeqFileParser(Options* options, const char* path, bool print, bool append = false, size_t BUFFER_SIZE = 4095);
	~SeqFileParser();

	inline int getFormat(){
		return _format;
	}

	inline void setFormat(int format){
		_format = format;
	}

	//get the next sequence from the file
	inline size_t getSeq(Sequence& seq, int numProcs, int myId, uint8_t prefixLength, int numBasesComp) {
		size_t ret;
		//int32_t numNs = 0, i;
		/*read the sequence from the file*/
		if (_format == FILE_FORMAT_FASTA) {
			if(_onlyOptical){
				Utils::exit("Distinguishing among optical and PCR sequences is only possible with FASTQ files\n");
			}
			ret = getFastaSeq(seq, numProcs, myId, prefixLength, numBasesComp);
		} else {
			ret = getFastqSeq(seq, numProcs, myId, prefixLength, numBasesComp);
		}
		return ret;
	}

	static inline void encode(uint8_t* s, size_t length) {
		uint8_t ch;
		for (size_t i = 0; i < length; i++) {
			ch = s[i];
			if (ch >= 'A' && ch <= 'Z') {
				ch -= 'A';
			} else if (ch >= 'a' && ch <= 'z') {
				ch -= 'a';
			} else {
				Utils::exit("Unexpected character %c at line %d in file %s\n",
						ch, __LINE__, __FILE__);
			}
			//s[i] = _codeTab[ch];
		}
	}

	inline void printSeq(Sequence *seq){
		if (_format == FILE_FORMAT_FASTA) {
			if(_onlyOptical){
				Utils::exit("Distinguishing among optical and PCR sequences is only possible with FASTQ files\n");
			}
			if(_compressFile){
				_printFastaSeqCompress(seq);
			} else {
				_printFastaSeqUncompress(seq);
			}
		} else {
			if(_compressFile){
				_printFastqSeqCompress(seq);
			} else {
				_printFastqSeqUncompress(seq);
			}
		}
	}
	void flushOutput();

private:
	/*private member functions*/
	void resizeBuffer(size_t nsize);
	int getFastaSeq(Sequence& seq, int numProcs, int myId, uint8_t prefixLength, int numBasesComp);
	int getFastqSeq(Sequence& seq, int numProcs, int myId, uint8_t prefixLength, int numBasesComp);

	void _printFastaSeqCompress(Sequence *seq);
	void _printFastqSeqCompress(Sequence *seq);
	void _printFastaSeqUncompress(Sequence *seq);
	void _printFastqSeqUncompress(Sequence *seq);

	/*buffered file operations*/
	inline int myfgetc(gzFile file) {
		/*check the end-of-file*/
		if (_fileBufferSentinel >= _fileBufferLength) {
			/*re-fill the buffer*/
			_fileBufferSentinel = 0;
			/*read file*/
			if(_compressFile){
				_fileBufferLength = gzread(file, _fileBuffer, 4096);
				if (_fileBufferLength == 0) {
					/*reach the end of the file*/
					if (gzeof(file)) {
						return -1;
					} else {
						Utils::exit("File reading failed in function %s line %d\n",
								__FUNCTION__, __LINE__);
					}
				}
			} else {
				_fileBufferLength = fread(_fileBuffer, 1, 4096, (FILE *)file);
				if (_fileBufferLength == 0) {
					/*reach the end of the file*/
					if (feof((FILE *)file)) {
						return -1;
					} else {
						Utils::exit("File reading failed in function %s line %d\n",
								__FUNCTION__, __LINE__);
					}
				}
			}
		}
		/*return the current character, and increase the sentinel position*/
		return _fileBuffer[_fileBufferSentinel++];
	}
	inline int myungetc(int ch, gzFile file) {
		if (_fileBufferSentinel >= 0) {
			_fileBuffer[--_fileBufferSentinel] = ch;
		} else {
			Utils::log("Two consecutive ungetc operations occurred\n");
			return -1; /*an error occurred, return end-of-file marker*/
		}
		return ch;
	}

private:
	bool _print;

	/*private member variables*/
	//buffer for file reading
	uint8_t* _buffer;
	size_t _length;
	size_t _size;

	//FASTA/FASTQ file handler
	int _compressFile;
	gzFile _fp;
	uint8_t* _fileBufferR;
	uint8_t* _fileBuffer;
	int _fileBufferLength;
	int _fileBufferSentinel;

	//
	int _format;
	bool _onlyOptical;

	static const uint8_t _codeTab[26];
	static const uint8_t _decodeTab[5];
	//static const uint8_t _complements[5];
};

#endif /* SEQFILEPARSER_H_ */
