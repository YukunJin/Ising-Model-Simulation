#include <vector>
#include <stack>
#include <list>
#include <map>
#include <iostream>
#include <random>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <thread>
using namespace std;
string decToBinary(int n)
{
    string ret = "";
    for (int i = 8; i >= 0; i--) {
        int k = n >> i;
        if (k & 1)
            ret += "1";
        else
            ret += "0";
    }
    return ret;
}

string arrayToBinary( vector<vector<int> > spins ) {
  std::stack<char> bin;
  string ret = "";
  for (int i=0 ; i<spins.size();i++){
    for (int j = 0; j<spins[i].size();j++){
    if(spins[i][j] == 1){
      bin.push('1');
    }
    else{
      bin.push('0');
    }
  }
}
  while(!bin.empty()){
    ret += bin.top();
    bin.pop();
  }
  return ret;
}
vector <vector<int > > randomFieldGenerator(int N){
  vector<vector<int >  > field;
  field.resize(N);
  int random;
  for (int i=0;i<N;i++){
    field[i].resize(N);
    for(int j=0;j<N;j++){
      random = rand()%2;
      if(random == 0){
        field[i][j] = -1;
      }
      else{
        field[i][j] = 1;
      }
    }
  }
  return field;
}
double calcEnergy(vector<vector<int > > &field,int Lx,int Ly, double J = 1 ){
  double energy  = 0.0;
  for (int i =0; i<Lx; i++){
    for(int j = 0 ; j<Ly ;j++){
      int S = field[i][j];
      int neighbor = field[(i+1)%Lx][j] + field[i][(j+1)%Ly]+
                    field[(Lx - i-1)%Lx][j] + field[i][(Ly-j-1)%Ly];
      energy += -S*neighbor;
    }
  }
  return energy*J;
}
double calcMag(vector<vector<int > > &field){
  double M = 0.0;
  int N = field.size();
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      int S = field[i][j];
      M += S;
    }
  }
  return pow((M/(N*N)),2);
}
double deltaE(vector<vector<int > > &field, int x, int y, double J = 1){
  int S = field[x][y];
  int N = field.size();
  int neighbor = field[(x+1)%N][y] + field[x][(y+1)%N]+
                field[(N-x-1)%N][y] + field[x][(N-y-1)%N];
  return 2*S*J*neighbor;
}

void mcmc( vector<vector<int > > &field , double beta_J){
  std::random_device rd;
  std::mt19937 mt(rd());
  int N = field.size();
  for (int i = 0; i < N*N ;i++){
    int x = rand()%N;
    int y = rand()%N;
    int s = field[x][y];
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(mt);
    double cost = deltaE(field,x,y);
    double probability = exp(-beta_J * cost);
    if(cost < 0){
      s *= -1;
    }
    else if (u <= probability){
      s *= -1;
    }
    field[x][y] = s;
  }
}
vector<vector<int > > coarse_grain(vector<vector<int > > &field ){
        vector<vector <int > > ret;
        vector<int> temp;
        int N = field.size();
        ret.resize(N/3);
        for (int i = 0;i<ret.size();i++){
            ret[i].resize(N/3);
            for (int j=0;j<ret[i].size();j++){
                ret[i][j] = 0;
                }
            }
    for(int m = 0 ; m<N/3;m++){
        for(int n = 0;n<N/3;n++){
            int up = 0;
            int down = 0;
            for (int i=m*3;i<m*3+3;i++){

                for (int j = n*3;j <n*3+3;j++){
                    if(field[i][j] == 1){
                        up++;
                        }
                    else{
                        down++;
                        }
                    }


                }
                if(up >= down){
                    ret[m][n] = 1;
                    }
                else{
                    ret[m][n] = -1;
                    }
            }
    }
            return ret;
    }
void main_func(vector<vector<int > > &field,double beta_J,string url){
  int N = field.size();
  double beta = beta_J;
  vector<double> EnergyList;
  std::vector<double> MagList;
  int count = 0;
  string folder  = to_string(N);
  ofstream file;
  string fileName = (url+"/beta_J=" + to_string(floorf(beta*100)/100) + ".txt");
  for (int sweep = 0; sweep < 10000; sweep++){
    count++;
    mcmc(field,beta);
    if(count > 100){
      double e_curr = calcEnergy(field,N,N);
      double m_curr = calcMag(field);
      EnergyList.push_back(e_curr);
      MagList.push_back(m_curr);
    }
  }
  file.open(fileName);
  while(file.is_open()){
    for(int i = 0; i < EnergyList.size();i++){
      file << "E "<< EnergyList[i];
      file << '\n';
    }

    for(int i = 0; i<MagList.size();i++){
      file << "M "<< MagList[i];
      file << '\n';
    }
    for(int i=0;i<field.size();i++){
      for (int j=0;j<field.size();j++){
        file << "S " << field[i][j];
        file << '\n';
      }
    }
    file.close();
}
}

void main_func_cg(std::vector<std::vector<int> > &field, double beta_J, string url){
  int C = field.size();
  vector<double> EnergyList_cg;
  std::vector<double> MagList_cg;
  ofstream file;
  string fileName = (url+"/beta_J=" + to_string(floorf(beta_J*100)/100) + ".txt");
  for (int sweep = 0; sweep<10000;sweep++){
      mcmc(field,beta_J);
      if(sweep > 100){
        std::vector<std::vector<int> > cg_field = coarse_grain(field);
        double e_curr = calcEnergy(cg_field,C/3,C/3);
        double m_curr = calcMag(cg_field);
        EnergyList_cg.push_back(e_curr);
        MagList_cg.push_back(m_curr);
      }
  }
  file.open(fileName);
  while(file.is_open()){
    for(int i = 0; i < EnergyList_cg.size();i++){
      file << "E "<< EnergyList_cg[i];
      file << '\n';
    }

    for(int i = 0; i<MagList_cg.size();i++){
      file << "M "<< MagList_cg[i];
      file << '\n';
    }
    for(int i=0;i<field.size();i++){
      for (int j=0;j<field.size();j++){
        file << "S " << field[i][j];
        file << '\n';
      }
    }
    file.close();
}

}

void writeToFile(std::vector<std::vector<int> > field,double beta_J,string url){
  ofstream file;

  string folder = to_string(field.size());
  string fileName = (url+"/"+folder+"/beta_J=" + to_string(floorf(beta_J*100)/100) + ".txt");
  file.open(fileName);
  for(int i=0;i<field.size();i++){
    for (int j=0;j<field.size();j++){
      file << "S " << field[i][j];
      file << '\n';
    }
  }
  file.close();
}


//////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
//////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
//////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
//////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////


int main(){
  int N = 27;
  int C = 81;
  vector<double> beta_J = {0.0,0.3,0.4,0.5,0.6,100};
  std::vector<std::vector<int> > native_field = randomFieldGenerator(N);
  std::vector<std::vector<int> > large_field;
  string url_native = "RENORMALIZED/native_field";
  string url_cg = "RENORMALIZED/cg_field";
  for (int i = 0;i<beta_J.size();i++){
     native_field = randomFieldGenerator(N);
     main_func(native_field,beta_J[i],url_native);
     //large_field = randomFieldGenerator(C);
    //thread th1(main_func,ref(native_field),beta_J[i],url_native);
    //thread th2(main_func_cg,ref(large_field),beta_J[i],url_cg);
    //th1.join();
    //th2.join();
    native_field.clear();
    //large_field.clear();
}
return 0;
}
