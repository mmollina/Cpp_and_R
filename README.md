# Very Basic Introduction to C++ anbd R

##Hello world

``` c++ 
#include <iostream>
int main()
{
    std::cout << "Hello, world!\n";
}
```

or yet

``` c++ 
#include <iostream>
using namespace std;
int main()
{
    cout << "Hello, world!\n";
}
```

Compiling

``` bash
g++ hello.cpp -O2 -Wformat -o hello
```

Another exemple using the `for` loop

``` c++ 
#include <iostream>
using namespace std;
int main()
{
  int a;
  int b;
  cout << "Enter the first number: ";
  cin >> a;
  cout << "\n";
  cout << "Enter the second number: ";
  cin >> b;
  cout << "\n";
  for(a; a < b; a++)
    {
        cout << a << "\n";
    }
}
```
Compiling

``` bash
g++ for_example.cpp -O2 -Wformat -o for_exemple
```

Another exemple: recombination fraction in a backcross


``` c++ 
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
  std::vector<std::vector<double> > rec(n_mar, std::vector<double>(n_mar, 0.0));
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

  for (int i = 0; i < n_mar; i++) 
    {
      for (int j = 0; j < n_mar; j++) 
	{
	cout << fixed;
	cout << setprecision(3);
	cout << rec[i][j] << " ";
      }
      cout << "\n";
    }
}
```
Compiling and running

``` bash
$ g++ bc_est.cpp -O2 -Wformat -o bc est
$ ./bcest
$ Enter the number of markers: 14
$ Enter the number of individuals: 103
$ Enter the name of the input file mouse.txt
```
R version

```{r}
dat<-read.table("mouse.txt")
id<-combn(1:ncol(dat),2)
rec<-apply(id, 2, function(x) sum(table(dat[,x[1]],dat[,x[2]])[2:3])/nrow(dat))
```

More realistc example
```{r}
source("simulate_diploid_populations.R")
require(onemap)
dat.bc<-sim.pop.bc(n.ind = 250, n.mrk = 500, ch.len = 200, missing = 0, n.ch = 1, verbose = FALSE)
dat.bc
dat<-dat.bc$geno
write.table(x=dat, file = "fake_bc.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

#using for
rec<-matrix(NA,ncol(dat), ncol(dat))
for(i in 1:(ncol(dat)-1)){
  for(j in (i+1):ncol(dat)){
    rec[i,j]<-sum(table(dat[,i], dat[,j])[2:3])/nrow(dat)
  }
}
require(fields)
image.plot(rec, col=rev(tim.colors()))

```

Using C++

``` bash
$ ./bcest
$ Enter the number of markers: 500
$ Enter the number of individuals: 250
$ Enter the name of the input file: fake_bc.txt
```

```{r}
y<-read.table(file = "rec_cpp.txt")
image.plot(as.matrix(y), col=rev(tim.colors()) )
```

##Conecting R with C++
