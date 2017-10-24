#ifndef KVPAR_H
#define KVPAR_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <algorithm>
#include <map>
using namespace std;

class kvpar
{
 public:
  kvpar(string f);
  string & getVal(string);
  void getVal(string, string &);
  void getVal(string, int &);
  void getVal(string, double &);
  void getVal(string, float &);
  void getVal(string, vector<string> &);
  void getVal(string, vector<int> &);
  void getVal(string, vector<double> &);
  void getVal(string, vector<float> &);
  bool getValAsBool(string);
  void getValAsRange(string, vector<int> &);

  void getVal(string, double *&);
  void getVal(string, int *&);
  void getVal(string, double *&, int &);
  void getVal(string, int *&, int &);

  void getFile(string, double *&, int&, int&);
  void getFile(string, double *&);

  void getFile(string, int *&, int&, int&);
  void getFile(string, int *&);

 private:
  map<string, string> kvMap;

};
#endif
