#ifndef GRAPHREADER_H
#define GRAPHREADER_H

using std::string;
using std::vector;
using std::size;
using std::endl;
using std::stringstream;
using std::istream;
using std::ios;
using std::to_string;

bool isLittleEndian(string s);

Graph parseGraph6String(string graphString);
Graph parsePregraphString(istream& pregraphString, bool littleEndian);
vector<Graph> parsePregraphs(string& s);
vector<Graph> parseMultiGraphs(vector<int>& ns);
vector<vector<int>> parseOrbit(string& s);


#endif //GRAPHREADER_H
