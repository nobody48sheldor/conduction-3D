#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include<string>

using namespace std;






// general functions
int write_data(vector<double> data, int n)
{
  ofstream output_file;
  output_file.open("data.txt");
  for (int i=0; i < n*n*n; i++)
  {
    output_file << data[i] << '/';
  }
  output_file.close();
  return 0;
}

auto linearized(vector<vector<vector<double>>> array, int n)
{
  vector<double> linear;
  linear.resize(n*n*n);

  for (int z=0; z<n; z++)
  {
    for (int y=0; y<n; y++)
    {
      for (int x=0; x<n; x++)
      {
        linear[z*n*n + y*n + x] = array[z][y][x];
      }
    }
  }
  return linear;
}

auto linspace(double l_0, double l_1, int n)
{
    vector<double> array;
    array.resize(n);

    for (int i = 0; i < n; i++)
    {
        array[i] = l_0 + (i*(l_1 - l_0)/n);
    }
    return array;
}







// making the axis (not necessary right now)

auto make_axis(int size, vector<double> L)
{
  vector<double> X;
  X.resize(size);

  vector<vector<double>> Y;
  Y.resize(size);

  vector<vector<vector<double>>> Z;
  Z.resize(size);

  for (int z=0; z<size; z++)
  {
    for (int y=0; y<size; y++)
    {
      for (int x=0; x<size; x++)
      {
        X[x] = L[x];
      }
      Y[y] = X;
    }
    Z[z] = Y;
  }
  return Z;
}

double grad(vector<vector<vector<double>>> tau, int x, int y, int z, double dl)
{
  double value = (tau[z][y][x+1] - 2*tau[z][y][x] + tau[z][y][x-1])/(dl*dl) + (tau[z][y+1][x] - 2*tau[z][y][x] + tau[z][y-1][x])/(dl*dl) + (tau[z+1][y][x] - 2*tau[z][y][x] + tau[z-1][y][x])/(dl*dl);
  return value;
}





// inital conditions

double tau_init_func(double x, double y, double z)
{
  double value = exp(- (x*x +y*y +z*z));
  return value;
}


auto tau_init(int size, vector<double> L)
{
  vector<vector<vector<double>>> array;
  array.resize(size);


  double value;

  for (int z=0; z<size; z++)
  {
    array[z].resize(size);
    for (int y=0; y<size; y++)
    {
      array[z][y].resize(size);
      for (int x=0; x<size; x++)
      {
        double X = L[x];
        double Y = L[y];
        double Z = L[z];
        array[z][y][x] = tau_init_func(X, Y, Z);
      }
    }
  }
  return array;
}





// solving function

auto solve(vector<vector<vector<double>>> tau_ti, int size, double alpha, double dl, double dt)
{
  vector<vector<vector<double>>> tau_t_next;
  tau_t_next.resize(size);
  double value;

  // making the empty array

  for (int z=0; z<size; z++)
  {
    tau_t_next[z].resize(size);
    for (int y=0; y<size; y++)
    {
      tau_t_next[z][y].resize(size);
    }
  }


  for (int z=1; z<size-1; z++)
  {
    cout << "z = " << z << endl;
    for (int y=1; y<size-1; y++)
    {
      for (int x=1; x<size-1; x++)
      {
        value = tau_ti[z][y][x] + alpha*(grad(tau_ti, x ,y, z, dl))*dt;
        tau_t_next[z][y][x] = value;
      }
      tau_t_next[z][y][size] = tau_t_next[z][y][size-1];
      tau_t_next[z][y][0] = tau_t_next[z][y][1];
    }
    tau_t_next[z][size-1] = tau_t_next[z][size-2];
    tau_t_next[z][0] = tau_t_next[z][1];
  }
  cout << "here " << endl;
  cout << size-1 << endl;
  tau_t_next[size-1] = tau_t_next[size-2];
  cout << "here " << endl;
  tau_t_next[0] = tau_t_next[1];
  return tau_t_next;
}





// main function

int main()
{
  // defining some values
  
  int n = 10;
  int nt = 5000;
  double limite = 5;
  double tmax = 35;
  double alpha = 1;


  // making some of the initial variable

  vector<double> L = linspace(-limite, limite, n);
  double dl = L[1] - L[0];
  vector<double> T = linspace(0, tmax, nt);
  double dt = T[1] - T[0];

  cout << "Making the initialisation ..." << endl << endl;
  cout << "   Calculating the initial values ..." << endl << endl;

  vector<vector<vector<double>>> tau_initial = tau_init(n,L);

  cout << "       --> Done !"<< endl << endl << endl;

  cout << "   fixing size of vectors <-- RAM expensive operation"<< endl << endl;

  vector<vector<vector<vector<double>>>> tau;
  tau.resize(nt);

  for (int t=0; t<4; t++)
  {
    tau[t].resize(n);
    for (int z=0; z<n; z++)
    {
      tau[t][z].resize(n);
      for (int y=0; y<n; y++)
      {
        tau[t][z][y].resize(n);
      }
    }
  }
  cout << "       --> Done !"<< endl << endl;

  cout << "initialisation done, ready to make some calculations" << endl << endl << endl;



  // initial condition put into the first frame

  tau[0] = tau_initial;


  // looping and making the next frame thanks to the previous one

  for (int t=1; t<nt; t++)
  {
    cout << "t = " << t << endl << endl;
    tau[t%4] = solve(tau[(t-1)%4], n, alpha, dl, dt);
    cout << "writting it up" << endl << endl;
    write_data(linearized(tau[t%4], n), n);
  }

  return 0;
}