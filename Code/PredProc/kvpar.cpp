#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
using namespace std;

#include "kvpar.h"

void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos)
    {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

kvpar::kvpar(string f)
{
  ifstream infile(f.c_str(), ios::in);
  if ( !infile ) {
    cerr << "Parameter file could not be opened." << endl;
    exit(1);
  }
  
  string temp, key, value;
  temp = "";


  while (getline(infile,temp)){
    istringstream input (temp);
    
    if (temp[0] != '#') {
      input >> key >> value;
      kvMap.insert(pair<string, string> (key,value));
    }
  }

}

string& kvpar::getVal(string k)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }
  return (*it).second;
}

void kvpar::getVal(string k, string & a)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }
  a = (*it).second;
}

void kvpar::getVal(string k, int & a)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }
  a = atoi(((*it).second).c_str());
}

void kvpar::getVal(string k, float & a)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }
  a = atof(((*it).second).c_str());
}

void kvpar::getVal(string k, double & a)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }
  a = atof(((*it).second).c_str());
}

void kvpar::getVal(string k, vector<string> & v)
{

  map<string, string>::iterator it;
  it = kvMap.find(k);

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  //string str = delSpaces((*it).second);
  //str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());

  tokenize((*it).second, v, ",");
}

void kvpar::getVal(string k, vector<int> & v)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  for(int i = 0; i < s.size(); i++)
    v.push_back(atoi(s[i].c_str()));
}

void kvpar::getValAsRange(string k, vector<int> & v)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ":");

  if(s.size() != 2 || atoi(s[1].c_str())<=atoi(s[0].c_str())){
    cout << "There is something wrong with the range associated with key: "  << k << endl;
    exit(1);
  }

  int i = atoi(s[0].c_str());
  while(i <= atoi(s[1].c_str())){
    v.push_back(i);
    i++;
  }

}

void kvpar::getVal(string k, vector<double> & v)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  for(int i = 0; i < s.size(); i++)
    v.push_back(atof(s[i].c_str()));
}

void kvpar::getVal(string k, vector<float> & v)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  for(int i = 0; i < s.size(); i++)
    v.push_back(atof(s[i].c_str()));
}



bool kvpar::getValAsBool(string k)
{

  map<string, string>::iterator it;
  it = kvMap.find(k);

  k = (*it).second;

  //Im sure there is a better way but for now...
  if(k =="yes"||k=="Yes"||k=="YES"||k=="t"||k=="T"||k=="true"||k=="True"||k=="TRUE"){
    return true;
  }
  else if(k =="no"||k=="No"||k=="NO"||k=="f"||k=="F"||k=="false"||k=="False"||k=="FALSE"){
    return false;
  }
  else{
    cout << "the true/false statement '" << k << "' is not recognized!" << endl;
    exit(1);
  }

}


void kvpar::getVal(string k, double *& v)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  v = new double[s.size()];

  for(int i = 0; i < s.size(); i++)
    v[i] = atof(s[i].c_str());

}


void kvpar::getVal(string k, int *& v)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  v = new int[s.size()];

  for(int i = 0; i < s.size(); i++)
    v[i] = atoi(s[i].c_str());

}




void kvpar::getVal(string k, double *& v, int & l)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  v = new double[s.size()];

  for(int i = 0; i < s.size(); i++)
    v[i] = atof(s[i].c_str());

  l = s.size();
}


void kvpar::getVal(string k, int *& v, int & l)
{
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  v = new int[s.size()];

  for(int i = 0; i < s.size(); i++)
    v[i] = atoi(s[i].c_str());

  l = s.size();
}



void kvpar::getFile(string k, double *& a, int & nrow, int & ncol){

  //get file name
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  nrow = atoi(s[1].c_str());
  ncol = atoi(s[2].c_str());

  cout << "Reading nrow=" << nrow << " by ncol=" << ncol << " file='"<< s[0] << "'" << endl;
 
  a = new double[nrow*ncol];

  ifstream f(s[0].c_str(), ios::in);
  if(!f){
    cout << "file '" << s[0] << "' cannot be opened" << endl;
    exit(1);
  }

  for(int j = 0; j < nrow; j++){
    for(int i = 0; i < ncol; i++){
      f >> a[i*nrow+j];  
    }
  }

  f.close();

}

void kvpar::getFile(string k, double *& a){

  //get file name
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  int nrow = atoi(s[1].c_str());
  int ncol = atoi(s[2].c_str());

  cout << "Reading nrow=" << nrow << " by ncol=" << ncol << " file='"<< s[0] << "'" << endl;

  a = new double[nrow*ncol];

  ifstream f(s[0].c_str(), ios::in);
  if(!f){
    cout << "file '" << s[0] << "' cannot be opened" << endl;
    exit(1);
  }

  for(int j = 0; j < nrow; j++){
    for(int i = 0; i < ncol; i++){
      f >> a[i*nrow+j];  
    }
  }

  f.close();

}



void kvpar::getFile(string k, int *& a, int & nrow, int & ncol){

  //get file name
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  nrow = atoi(s[1].c_str());
  ncol = atoi(s[2].c_str());

  cout << "Reading nrow=" << nrow << " by ncol=" << ncol << " file='"<< s[0] << "'" << endl;
 
  a = new int[nrow*ncol];

  ifstream f(s[0].c_str(), ios::in);
  if(!f){
    cout << "file '" << s[0] << "' cannot be opened" << endl;
    exit(1);
  }

  for(int j = 0; j < nrow; j++){
    for(int i = 0; i < ncol; i++){
      f >> a[i*nrow+j];  
    }
  }

  f.close();

}

void kvpar::getFile(string k, int *& a){

  //get file name
  map<string, string>::iterator it;
  it = kvMap.find(k);
  vector<string> s;

  if(it == kvMap.end()){
    cout << "Key-value pair not found for key: " << k << endl;
    exit(1);
  }

  tokenize((*it).second, s, ",");

  int nrow = atoi(s[1].c_str());
  int ncol = atoi(s[2].c_str());

  cout << "Reading nrow=" << nrow << " by ncol=" << ncol << " file='"<< s[0] << "'" << endl;

  a = new int[nrow*ncol];

  ifstream f(s[0].c_str(), ios::in);
  if(!f){
    cout << "file '" << s[0] << "' cannot be opened" << endl;
    exit(1);
  }

  for(int j = 0; j < nrow; j++){
    for(int i = 0; i < ncol; i++){
      f >> a[i*nrow+j];  
    }
  }

  f.close();

}
