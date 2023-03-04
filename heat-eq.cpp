#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include<string>

using namespace std;






// general functions
int write_data(vector<double> data, int n, int index)
{
  ofstream output_file;
  output_file.open("data/data_" + to_string(index) + ".txt");
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

double average(vector<vector<vector<double>>> array, int size)
{
  double avg = 0;

  for (int z=0; z<size; z++)
  {
    for (int y=0; y<size; y++)
    {
      for (int x=0; x<size; x++)
      {
        avg += array[z][y][x];
      }
    }
  } 
  return avg/(size*size*size);
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


// reduce size of array

auto red_size(vector<vector<vector<double>>> array, int n)
{
  vector<vector<vector<double>>> arr_red;
  //resizing
  
  arr_red.resize(n/2);
  for (int z=0; z<(n/2); z++)
  {
    arr_red[z].resize(n/2);
    for (int y=0; y<(n/2); y++)
    {
      arr_red[z][y].resize(n/2);
    }
  }

  // making the average of 4 point of array to reduce its size


  
  for (int z=0; z<(n/2); z++)
  {
    for (int y=0; y<(n/2); y++)
    {
      for (int x=0; x<(n/2); x++)
      {
        arr_red[z][y][x] = (1.0/8.0) * (array[2*z][2*y][2*x] + array[2*z][2*y][2*x + 1] + array[2*z][2*y + 1][2*x + 1] + array[2*z][2*y + 1][2*x] + array[2*z + 1][2*y][2*x] + array[2*z + 1][2*y][2*x + 1] + array[2*z + 1][2*y + 1][2*x + 1] + array[2*z + 1][2*y + 1][2*x]);
      }
    }
  }

  return arr_red;
}



// inital conditions

// function for the initial heat distribution

double tau_init_func(double x, double y, double z)
{
  double value;
  // value = 20*exp(- 0.2*(x*x +y*y +z*z));
  if (x < 0)
  {
    value = 30+(-7/200)*x;
  }
  else
  {
    if (y*y + z*z < 100*100)
    {
      if (x < 100)
      {
        value = -7;
      }
      else
      {
        value = 20;
      }
    }
    else
    {
      value = 20;
    }
  }
  return value;
}


// function for the alpha coefficient

double alpha_func(double x, double y, double z)
{
  double value;
  // value = 20*exp(- 0.2*(x*x +y*y +z*z));
  if (x < 0)
  {
    value = 12;
  }
  else
  {
    if (y*y + z*z < 5*5)
    {
      if (x > 5)
      {
        value = 111;
      }
      else
      {
        value = 1.5;
      }
    }
    else
    {
      value = 18.46;
    }
  }
  return value;
}


// making a 3D array with a "func" function

auto array_3D(int size, vector<double> L, double (*func)(double, double, double))
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
        array[z][y][x] = func(X, Y, Z);
      }
    }
  }
  return array;
}





// solving function

auto solve(vector<vector<vector<double>>> tau_ti, int size, vector<vector<vector<double>>> alpha, double dl, double dt)
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


  // solving with numerical methods

  for (int z=1; z<size-1; z++)
  {
    for (int y=1; y<size-1; y++)
    {
      for (int x=1; x<size-1; x++)
      {
        value = tau_ti[z][y][x] + alpha[z][y][x]*(grad(tau_ti, x ,y, z, dl))*dt;
        tau_t_next[z][y][x] = value;
      }
      tau_t_next[z][y][size-1] = tau_t_next[z][y][size-2];
      tau_t_next[z][y][0] = tau_t_next[z][y][1];
    }
    tau_t_next[z][size-1] = tau_t_next[z][size-2];
    tau_t_next[z][0] = tau_t_next[z][1];
  }
  tau_t_next[size-1] = tau_t_next[size-2];
  tau_t_next[0] = tau_t_next[1];
  return tau_t_next;
}





// main function

int main()
{
  // defining some values
  
  int n = 4*12;
  int nt = 400;
  double limite = 200;
  double tmax = 50;


  // making some of the initial variable

  vector<double> L = linspace(-limite, limite, n);
  double dl = L[1] - L[0];
  vector<double> T = linspace(0, tmax, nt);
  double dt = T[1] - T[0];


  vector<vector<vector<double>>> alpha = array_3D(n,L,alpha_func);

  cout << "Making the initialisation ..." << endl << endl;
  cout << "   Calculating the initial values ..." << endl << endl;

  vector<vector<vector<double>>> tau_initial = array_3D(n,L,tau_init_func);
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

  vector<vector<vector<double>>> reduced;
  reduced.resize(n/2);
  for (int z=0; z<n/2; z++)
  {
    reduced[z].resize(n/2);
    for (int y=0; y<n/2; y++)
    {
      reduced[z][y].resize(n/2);
    }
  }


  tau[0] = tau_initial;
  double avg = average(tau[0],n);
  cout << "average = " << avg << endl;

  reduced = red_size(tau[0],n);
  write_data(linearized(red_size(reduced,n/2),n/4),n/4, 0);
  system("python3 plot-a-frame.py 0");


  // looping and making the next frame thanks to the previous one

  for (int t=1; t<nt; t++)
  {
    cout << "t = " << t << endl;
    tau[t%4] = solve(tau[(t-1)%4], n, alpha, dl, dt);
    double avg = average(tau[t%4],n);
    cout << "average = " << avg << endl;
    reduced = red_size(tau[t%4], n);
    write_data(linearized(red_size(reduced,n/2), n/4), n/4, t);
    string cmd;
    cmd = "python3 plot-a-frame.py " + to_string(t) + " &";
    const char *command = cmd.c_str();
    system(command);
  }

  return 0;
}