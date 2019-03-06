#include <vector>
#include <fstream>
#include <iostream>
using namespace std;
class Ising{
  vector<vector<int > > spin;
  int Lx;
  int Ly;
  double beta;
  double J;
  double h;
public:
  Ising();
  Ising(&vector<vector<int > > NewSpin, newLx, newLy, newBeta, newJ, newH);
  vector<vector<int > > readSpinFile(string fileName);
  double calcEnergy();
  double deltaE(int spin_to_flip);
};

Ising(&vector<vector<int> > NewSpin, newLx, newLy, newBeta, newJ, newH){
  spin.resize(NewSpin.size());
  for (int i=0;i<spin.size();i++){
    spin[i].resize(NewSpin[i].size())
    for (int j=0;j<spin[i].size();j++){
      spin[i][j]  = NewSpin[i][j];
    }
  }
  Lx = newLx;
  Ly = newLy;
  beta = newBeta;
  J = newJ;
  h = newH;
}

vector<vector<int> > readSpinFile(string fileName){
  string line;
  ifstream spinsFile;
  spinsFile.open(filename);
  vector<vector<int> > spin;

  if(spinsFile.is_open()){
    while(getline(spinsFile,line)){
      spinsFile >> node_i >> h_value >> spin;
      setH(node_i,h_value);
      setSpin(node_i,spin);
    }
    spinsFile.close();
  }
  return true; 
}
