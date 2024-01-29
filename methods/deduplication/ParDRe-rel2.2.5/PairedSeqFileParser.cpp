#include "PairedSeqFileParser.h"


PairedSeqFileParser::PairedSeqFileParser(Options* options, const char* path1, const char* path2, size_t BUFFER_SIZE) {
	_parser1 = new SeqFileParser(options, path1, false, false, BUFFER_SIZE);
	_parser2 = new SeqFileParser(options, path2, false, false, BUFFER_SIZE);
}

PairedSeqFileParser::~PairedSeqFileParser() {
	delete _parser1;
	delete _parser2;
}
