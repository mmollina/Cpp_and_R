#include <iostream>
#include <vector>
#include <fstream> 
#include <string>
#include<iomanip>

using namespace std;
int main()
{


  int n_mar;
  int n_ind;
  double ct;
  cout << "Enter the number of markers: ";
  cin >> n_mar;
  cout << "\n";
  cout << "Enter the number of individuals: ";
  cin >> n_ind;
  cout << "\n";
  string filename;
  ifstream in;
  cout << "Enter the name of the input file ";
  cin >> filename;
  in.open( filename.c_str() );
  if (!in) {
    cout << "Cannot open file.\n";
    return (-1);
  }
  std::vector<std::vector<int> > geno(n_ind, std::vector<int>(n_mar, 0));
  std::vector<std::vector<double> > rec(n_mar, std::vector<double>(n_mar, 0.5));
  for (int i = 0; i < geno.size(); i++) {
    for (int j = 0; j < geno[1].size(); j++) {
      in >> geno[i][j];
    }
  }
  in.close();

  for (int i = 0; i < (n_mar-1); i++) 
    {
      for (int j = (i+1); j < n_mar; j++) 
	{
	  ct=0.0;
	  for(int k=0; k < n_ind; k++)
	    {
	      if(geno[k][i]!=geno[k][j])
		ct++;
	    }
	  rec[i][j]=ct/n_ind;
	}
    }

  ofstream fout("rec_cpp.txt"); //opening an output stream fo
  for (int i = 0; i < n_mar; i++) 
    {
      for (int j = 0; j < n_mar; j++) 
	{
	fout << rec[i][j] << " ";
      }
      fout << "\n";
    }
}
