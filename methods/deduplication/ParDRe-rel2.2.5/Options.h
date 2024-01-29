/*
 * Options.h
 *
 *  Created on: Jun 1, 2015
 *      Author: gonzales
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_
#include "Macros.h"
#include "Utils.h"
#include "mpi.h"

class Options
{
public:
	Options();
	virtual ~Options();

	/*virtual functions*/
	void printUsage();
	bool parse(int argc, char* argv[]);

	//get the input file list
	inline vector<string> getInputFileList() {
		return _readFileNames;
	}

	inline vector<string> getOutFileList() {
		return _outFileNames;
	}

	inline int getNumThreads() {
		return _numThreads;
	}
	inline bool isPaired() {
		return _paired;
	}
	inline uint8_t getMissMatch(){
		return _missMatch;
	}
	inline uint8_t getPrefixLen(){
		return _prefixLen;
	}
	inline int getNumBasesCompared(){
		return _numBasesCompared;
	}
	inline bool getPrint(){
		return _print;
	}
	inline uint32_t getBlockSize(){
		return _blockSize;
	}
	inline bool onlyOptical(){
		return (_opticalDist >= 0);
	}
	inline int getOpticalDist(){
		return _opticalDist;
	}
	inline int getCompressFile(){
		return _compressedFiles;
	}

private:
	void _setDefaults();
    string _makeOutputName(const string& inName);

	/*member variables*/
	int _compressedFiles;
	vector<string> _readFileNames;
	vector<string> _outFileNames;
	int _numThreads;
	bool _paired;
	bool _print;

	uint8_t _missMatch;
	uint8_t _prefixLen;
	int _numBasesCompared;

	uint32_t _blockSize;

	int _opticalDist;
};

#endif /* OPTIONS_H_ */
