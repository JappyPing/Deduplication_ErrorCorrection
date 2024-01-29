#include "SeqFileParser.h"
#include "Utils.h"
#include <zlib.h>

/*const uint8_t SeqFileParser::_codeTab[26] = { 0, 4, 1, 4, 4, 4, 2, //A -> G
		4, 4, 4, 4, 4, 4, 4, //H->N
		4, 4, 4, 4, 4, 3, 4, //O->U
		4, 4, 4, 4, 4 //V->Z
		};*/


SeqFileParser::SeqFileParser(Options* options, const char* path, bool print, bool append, size_t BUFFER_SIZE) {

	_print = print;

	_onlyOptical = options->onlyOptical();
	_compressFile = options->getCompressFile();

	if(!print){
		/*create the file buffer*/
		_fileBufferSentinel = 0;
		_fileBufferLength = 0;
		_fileBufferR = new uint8_t[4096 + 8];
		if (_fileBufferR == NULL) {
			Utils::exit("Memory allocation failed in file %s in line %d\n",
					__FUNCTION__, __LINE__);
		}
		_fileBuffer = _fileBufferR + 8; /*make it aligned*/

		/*open the input file*/
		if(_compressFile){
			_fp = gzopen(path, "rb");
		} else {
			_fp = (gzFile) fopen(path, "rb");
		}
		if (_fp == NULL) {
			Utils::cmpExit("Failed to open file: %s\n", path);
		}

		//allocate buffer for file reading
		_size = BUFFER_SIZE;
		_length = 0;
		_buffer = new uint8_t[_size + 1];
		if (_buffer == NULL) {
			Utils::exit("Memory allocation failed in file %s in line %d\n",
					__FUNCTION__, __LINE__);
		}

		//detecting the file format in the first line
		int ch;

		while ((ch = myfgetc(_fp)) != -1 && ch != '>' && ch != '@' && ch != '\n')
				;
			if (ch == -1 || ch == '\n') {
				Utils::exit("Unrecognized file format\n");
			} else if (ch == '>') {
				_format = FILE_FORMAT_FASTA;
				myungetc(ch, _fp);
				Utils::log("FASTA format identified\n");
			} else {
				_format = FILE_FORMAT_FASTQ;
				myungetc(ch, _fp);
				Utils::log("FASTQ format identified\n");
			}
	} else if (!options->getPrint()) {
		// no need to open output file, not printing output
		_fp = NULL;
	} else{
		/*open the output file*/
		if(_compressFile){
			if(append){
				_fp = gzopen(path, "a");
			} else {
				_fp = gzopen(path, "w");
			}
			if (_compressFile > 1) {
				gzsetparams(_fp, _compressFile, Z_DEFAULT_STRATEGY);
			}
		} else {
			if(append){
				_fp = (gzFile) fopen(path, "a");
			} else {
				_fp = (gzFile) fopen(path, "w");
			}
		}
		if (_fp == NULL) {
			Utils::exit("Failed to open output file: %s\n", path);
		}
	}
}

SeqFileParser::~SeqFileParser() {
	/*check the file format*/
	//close the file
	if (!_fp) {
		// no need to close, never opened
	}
	else if(_compressFile){
		gzclose(_fp);
	} else {
		fclose((FILE *)_fp);
	}

	//release the buffer
	if(!_print){
		if(_fileBufferR){
			delete[] _fileBufferR;
		}
		if(_buffer){
			delete[] _buffer;
		}
	}
}

void SeqFileParser::resizeBuffer(size_t nsize) {
	if (nsize <= _size) {
		return;
	}

	//allocate a new buffer
	_size = nsize * 2;
	uint8_t* nbuffer = new uint8_t[_size];
	if (!nbuffer) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	//copy the old data
	memcpy(nbuffer, _buffer, _length);

	//release the old buffer
	delete[] _buffer;
	_buffer = nbuffer;
}

int SeqFileParser::getFastaSeq(Sequence& seq, int numProcs, int myId, uint8_t prefixLength, int numBasesComp) {
	int ch;

	//find the header
	while ((ch = myfgetc(_fp)) != -1 && ch != '>')
		;
	if (ch == -1) {
		return 0; //reach the end of file
	}
	//read the sequence name (only one line)
	_length = 0;
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = ch;
	}
	if (ch == -1) {
		Utils::exit("Incomplete file\n");
	}
	_buffer[_length] = '\0';

	/*trim characters /[12]$ like BWA*/
	if (_length > 2 && _buffer[_length - 2] == '/'
			&& (_buffer[_length - 1] == '1' || _buffer[_length - 1] == '2')) {
		_length -= 2;
		_buffer[_length] = '\0';
	}

	//save the sequence name
	seq.setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) seq._name, (char*) _buffer);

	//read the sequence bases
	_length = 0;
	do {
		//filter out the blank lines
		while ((ch = myfgetc(_fp)) != -1 && (ch == '\r' || ch == '\n'))
			;
		if (ch == -1)
			break; //reach the end of file
		if (ch == '>') { //reaching another sequence
			myungetc(ch, _fp);
			break;
		}

		//encode and save the base
		if (ch >= 'A' && ch <= 'Z') {
			ch -= 'A';
		} else if (ch >= 'a' && ch <= 'z') {
			ch -= 'a';
		} else {
			Utils::exit("Unexpected character %c at line %d in file %d\n", ch,
					__LINE__, __FILE__);
		}
		//ch = _codeTab[ch];

		//save the current encoded base
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = ch;
	} while (1);

	// save the sequence prefix
	uint64_t key = seq.setKey(_buffer, prefixLength, _length);
	if((int) (key%numProcs) != myId){
		return -1;
	}

	// save the suffix
	seq.setSuffix(_buffer, numBasesComp);

	return _length;
}

int SeqFileParser::getFastqSeq(Sequence& seq, int numProcs, int myId, uint8_t prefixLength, int numBasesComp) {
	int ch;

	//find the header
	while ((ch = myfgetc(_fp)) != -1 && ch != '@')
		;
	if (ch == -1)
		return 0; //reach the end of file

	//read the sequence name (only one line)
	_length = 0;
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = ch;
	}
	if (ch == -1) {
		Utils::exit("Incomplete file\n");
	}
	_buffer[_length] = '\0';

	/*trim characters /[12]$ like BWA*/
	/*if (_length > 2 && _buffer[_length - 2] == '/'
			&& (_buffer[_length - 1] == '1' || _buffer[_length - 1] == '2')) {
		_length -= 2;
		_buffer[_length] = '\0';
	}*/

	//save the sequence name
	seq.setNameSize(_length + 1); /*adjust the name buffer size*/
	strcpy((char*) seq._name, (char*) _buffer);

	// Get the coordinates if only interested in optical reads
	if(_onlyOptical){
		seq.setOpticalCoords();
	}

	//read the sequence bases
	_length = 0;
	do {
		//filter out the blank lines
		while ((ch = myfgetc(_fp)) != -1 && (ch == '\r' || ch == '\n'))
			;
		if (ch == -1)
			Utils::exit("Incomplete FASTQ file\n");
		if (ch == '+')
			break; //the comment line

		//encode and save the base
		/*base space encoding*/
		if (ch >= 'A' && ch <= 'Z') {
			ch -= 'A';
		} else if (ch >= 'a' && ch <= 'z') {
			ch -= 'a';
		} else {
			Utils::exit("Unexpected character %c at line %d in file %s\n",
						ch, __LINE__, __FILE__);
		}
		//ch = _codeTab[ch];

		//save the current encoded base
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		_buffer[_length++] = ch;
	} while (1);

	if(!_length){
		return 0;
	}

	// save the sequence prefix
	uint64_t key = seq.setKey(_buffer, prefixLength, _length);
	if((int) (key%numProcs) == myId){
		// save the suffix
		seq.setSuffix(_buffer, numBasesComp);
	}

	//read the comment line (only one line)
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n')
		;
	if (ch == -1)
		Utils::exit("Incomplete FASTQ file\n");

	//read the quality scores
	_length = 0;
	while ((ch = myfgetc(_fp)) != -1 && ch != '\n') {
		if (_length >= _size) {
			resizeBuffer(_size + 256);
		}
		if (ch >= 33 && ch <= 127) {
			_buffer[_length++] = ch;
		}

		if (_length > seq._length)
			break;
	}

	if(!_length){
		return 0;
	}

	/*for base-space reads*/
	if (seq._length != _length) {
		Utils::exit(
				"The number of bases is not equal to the number of quality scores\n");
	}

	if((int) (key%numProcs) != myId){
		return -1;
	}

	seq.setQualsSize(_length);
	memcpy(seq._quals, _buffer, _length);

	return seq._length;
}

void SeqFileParser::_printFastaSeqUncompress(Sequence *seq){
	fputc('>', (FILE *)_fp);
	seq->printUncompress((FILE *)_fp, false);
}

void SeqFileParser::_printFastaSeqCompress(Sequence *seq){
	gzputc(_fp, '>');
	seq->printCompress(_fp, false);
}

void SeqFileParser::_printFastqSeqUncompress(Sequence *seq){
	fputc('@', (FILE *) _fp);
	seq->printUncompress((FILE *)_fp, true);
}

void SeqFileParser::_printFastqSeqCompress(Sequence *seq){
	gzputc(_fp, '@');
	seq->printCompress(_fp, true);
	//gzflush(_fp, Z_FINISH);
}

void SeqFileParser::flushOutput() {
	if (_compressFile) {
		gzflush(_fp, Z_FINISH);
	}
	else {
		fflush((FILE *)_fp);
	}
}
