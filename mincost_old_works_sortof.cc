//Author: Mark Hubenthal
//Date: 12/08/2008
//Program: mincost.cc
//Solves the min cost flow problem given integer arc capacities, unit costs, and divergence values using the epsilon relaxation algorithm with epsilon scaling
//To compile: g++ mincost.cc -o mincost
//To run: ./mincost file.dat
//Make sure that file.dat is formatted correctly

#include <iostream>

#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <deque>
#include <time.h>

using namespace std;

double minpot1( vector<int> v, double *r ){
  double m;
  m = -r[v.back()];
  v.pop_back();

  while(v.empty() != true){
    m = min(m, -r[v.back()]);
    v.pop_back();
  }
  return m;
}

double minpot2( vector<int> v, double *r ){
  double m;
  
  m = r[v.back()];
  v.pop_back();

  while(v.empty() != true){
    m = min(m, r[v.back()]);
    v.pop_back();
  }
  return m;
}

//initializes x given a potential u. r is related to u
//by the equation r_j = d_j + u_i - u_i' where j = (i, i') is an arc
void initflow( int *x, double *r, int A, int *c){
  int i;
  for(i=1; i<=A; i++){
    if( r[i] >= 0){
      x[i] = 0;
    }else{
      x[i] = c[i];
    }
  }
}

//computes the cost of a given flow
double cost( int *x, int *d, int A){
  int i;
  double totalcost = 0;
  for(i = 1; i<=A; i++){
    totalcost += double(x[i])*double(d[i]);
  }
  return totalcost;
}

//checks whether a given flow x is feasible with respect to capacities
//and whether s = 0
bool checkfeas( int *x, int *c, int *s, int A, int N){
  int i;
  bool feasible=true;
  
  for(i=1; i<=A; i++){
    if( (x[i] > c[i]) || (x[i]<0) ){ feasible = false; }
  }
  for(i=1; i<=N; i++){
    if(s[i] != 0){ feasible = false; }
  }
  return feasible;
}

void printflow( int *x, int A){
  int i;
  for(i=1; i<=A; i++){
    cout<<"x["<<i<<"] = "<<x[i]<<endl;
  }
}

int main( int argc, char *argv[] ) {

//   ifstream inputdata( argv[1] );
  ifstream inputdata;
  int databuffer; //variable storing the last read integer from the file.dat
  int i; int j;
  time_t t_i; //start time of program
  time_t t_f; //end time of program
  double dif; //run time of program

  time(&t_i); //records start time of program

  //Reading in the number of nodes and arcs from the first line of file.dat
  inputdata.open(argv[1]);
  inputdata >> databuffer;
  const int N = databuffer;
  inputdata >> databuffer;
  const int A = databuffer;

  int startnode[A];
  int endnode[A];
  int d[A];
  double maxd; //maximum of absolute values elements in d
  int c[A]; //upper capacity on each arc
  int b[N]; //supply at each node
  int degree[N]; //degree of each node
  int maxdeg;
  int maxdegnode;
  for(j=1; j <= N; j++){ degree[j] = 0; } //initialize degree of each node to zero

  //properly import the data into the arrays startnode[], endnode[], c[], d[]
  
  i=1;
  while(i <= A ){
    inputdata >> databuffer;
    startnode[i] = databuffer;
    degree[databuffer]++;

    inputdata >> databuffer;
    endnode[i] = databuffer;
    degree[databuffer]++;

    inputdata >> databuffer;
    d[i] = databuffer;
    inputdata >> databuffer;
    c[i] = databuffer;
    i++;
  }

  //import supply data into b[]
  i = 1;
  while(i <= N ){
    inputdata >> databuffer;
    b[i] = databuffer;
    i++;
    
  }
  inputdata.close();

  //Compute the maximum degree of any given node
  maxdeg=0;
  for(j=1; j<=N; j++){ 
    if(maxdeg < degree[j]){ 
      maxdegnode = j;
      maxdeg = degree[j];
    }
  }



  int arcout[N][maxdeg]; //Records the indices from the array startnode[] corresponding to arcs out of each node

  int arcin[N][maxdeg]; //Records the indices from the array endnode[] corresponding to arcs into each node

  int numarcout[N];
  int numarcin[N];

  //initialize so that there are zero arcs coming in or out of each node
  for(i=1; i <= N; i++){
    numarcout[i] = 0;
    numarcin[i] = 0;
  }

  //build arcout and arcin
  for(i=1; i <= A; i++){
    j=numarcin[endnode[i]];
    numarcin[endnode[i]] = j+1;
    arcin[endnode[i]][numarcin[endnode[i]]] = i;
    j = numarcout[startnode[i]];
    numarcout[startnode[i]] = j+1;
    arcout[startnode[i]][numarcout[startnode[i]]] = i;
  }


  //compute the maximum of |d_{j}|
  maxd = fabs(double(d[1]));
  for(i=1; i<=A; i++){
    if(fabs(double(d[i]))>maxd){ 
      maxd = fabs(double(d[i])); 
    }
  }

  double e = maxd;
  double u[N];
  int x[A];
  double r[A]; //reduced cost d_j + u_i - u_i' for each arc
  int s[N]; //excess supply at each node, s = Ex - b
  deque<int> nodeq;
  int ibar;

  //Start with zero node potentials
  for(i=1; i<=N; i++){
    u[i] = 0;
  }
  //Choose x to be feasible and initialize reduced cost r
  for(i=1; i<=A; i++){
    r[i] = double(d[i]) + u[startnode[i]] - u[endnode[i]];
    if( r[i] >= 0){
      x[i] = 0;
    }else{
      x[i] = c[i];
    }
  }

  //Compute s
  for(i=1; i<=N; i++){

    s[i] = -b[i];
    for(j=1; j<=numarcout[i]; j++){
      s[i] = s[i] + x[arcout[i][j]];
    }
    for(j=1; j<=numarcin[i]; j++){
      s[i] = s[i] - x[arcin[i][j]];
    }
  }

  //Put all nodes with positive s in the queue
  for(i=1; i<=N; i++){
    if(s[i] > 0){ 
      nodeq.push_back(i);
    }
  }

  //---------------------------------------------------------------------//
  //---------------------------------------------------------------------//
 
  //Main part of algorithm

  bool proceed;
  int k;
  double alpha;
  double alpha1;
  double alpha2;
  int beta;
  vector<int> temp;
  int arc;

  while(e >= 1/( double(N) )){  
    
    //change flow so it satisfies epsilon complemetary slackness with respect
    //to new epsilon.
    for(i=1; i<=A; i++){
      if( r[i] >= 0){
    x[i] = 0;
      }else{
    x[i] = c[i];
      }
    }

  //   Compute s
    for(i=1; i<=N; i++){
      s[i] = -b[i];
    for(j=1; j<=numarcout[i]; j++){
      k = arcout[i][j];
      s[i] += x[k];
    }
    for(j=1; j<=numarcin[i]; j++){
      k = arcin[i][j];
      s[i] -= x[k];
    } 
    }
    
    //Put all nodes with positive s in queue
    for(i=1; i<=N; i++){
      if(s[i] > 0){ 
    nodeq.push_back(i);
      }
    }

    while( nodeq.empty() == false){
  
      ibar = nodeq.back();
      nodeq.pop_back();
    
      while(s[ibar] > 0){

    proceed = true;
    temp.clear();

    for(j=1; j<=numarcout[ibar] && proceed == true; j++){
      k = arcout[ibar][j];

      if( (r[k] <= e) && (r[k] >= (e/2)) && (x[k] > 0) ){ 
        arc = k; 
        proceed = false;
      }
    }

    if(proceed == false ){
      //if proceed == false, then there are black arcs out of ibar
      k = arc;
      beta = min( x[k], s[ibar] );
      x[k] = x[k] - beta;
      s[ibar] = s[ibar] - beta;
      s[endnode[k]] = s[endnode[k]] + beta;

      if( (s[endnode[k]] > 0) && (s[endnode[k]] - beta <= 0) ){ 
        nodeq.push_front(endnode[k]);
      }

    }else{
      //This portion executes if there are no black arcs out of ibar
      for( j=1; j<=numarcin[ibar] && proceed==true; j++ ){
      
        k = arcin[ibar][j];

        if( (r[k] >= (-e)) && (r[k] <= (-e/2)) && (x[k] < c[k]) ){
          arc = k;
          proceed = false;
          //Checks for white arcs into ibar. 
          //Loop stops after first arc is found
        }
      }

      if(proceed == false ){
        k = arc;
        beta = min( c[k] - x[k], s[ibar]);
        x[k] = x[k] + beta;
        s[ibar] = s[ibar] - beta;
        s[startnode[k]] = s[startnode[k]] + beta;

        if( (s[startnode[k]] > 0) && (s[startnode[k]] - beta) <= 0 ){
          nodeq.push_front(startnode[k]);
        }

      }
    }

    //Increase potential u[ibar] if all arcs are red

    if( proceed==true ){  
      //the vector temp will first store the labels of all arcs out of
      //ibar with positive flow
      for(j=1; j <= numarcout[ibar]; j++){
        k = arcout[ibar][j];
        if( x[k] > 0 ){ temp.push_back(k); }
      }

      if(temp.empty() == false){
        alpha1 = minpot1( temp, r );
        temp.clear();
        //minpot1(temp, r) considers all arcs in temp, which are arcs
        //out of ibar with positive flow. Then it finds the minimum
        //of -r[j] over all such arcs. Finally we clear the vector temp
        //for future use.

        for(j=1; j <= numarcin[ibar]; j++){
          //Now we use temp to store the labels of all arcs going
          //into ibar with flow strictly less than the upper capacity.
          k = arcin[ibar][j];
          if( x[k] < c[k] ){ temp.push_back(k); }
        }
        if(temp.empty() == false){
          alpha2 = minpot2( temp, r );
          temp.clear();
          alpha = min(alpha1,alpha2) + e;
          //minpot2(temp, r) minimizes r[j] over all arcs j in temp,
          //which are arcs into ibar with flow strictly less than
          //the upper capacity.
          //We then take alpha to be the minimum of alpha1 and alpha2 + e.

        }else{
          alpha = alpha1 + e;
          //in the case that there are no arcs into ibar with
          //flow strictly less than the upper capacity.
        }
      }else{
        for(j=1; j <= numarcin[ibar]; j++){
          k = arcin[ibar][j];
          if( x[k] < c[k] ){ temp.push_back(k); }
        }
        alpha = minpot2( temp, r) + e;
        temp.clear();
      }

      //update u
      u[ibar] = u[ibar] + alpha;

      //update r
      for(j=numarcout[ibar]; j >= 1; j--){
        k = arcout[ibar][j];
        r[k] = r[k] + alpha;
      }
      for(j=numarcin[ibar]; j >= 1; j--){
        k = arcin[ibar][j];
        r[k] = r[k] - alpha;
      }
    }     
      }
    }

    //scale epsilon down by 1/2
    e = e/2;
  }
  
  cout.setf(ios::floatfield);
  cout.precision(20); //sets number of digits to print out
  
  time(&t_f); //records end time of program
  dif = difftime(t_f, t_i);
  cout<<"Run time: "<<dif<<" seconds"<<endl;

  //return cost associated with obtained flow x
  if( checkfeas(x,c,s,A,N) ){ 
    cout<<"Current flow is feasible"<<endl;
    cout<<"Minimum cost is: "<<cost(x,d,A)<<endl;
  }else{
    cout<<"Current flow not feasible"<<endl;
  }

}
        
      
       
