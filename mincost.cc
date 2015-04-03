//Author: Mark Hubenthal
//Date: 12/2008 (revised 2015)
//Program: mincost.cc
//
//Solves the min cost flow problem given integer arc capacities, unit costs, 
//and divergence values using the epsilon relaxation algorithm with epsilon 
//scaling. In simple terms, it finds a flow vector
//x indexed by the arcs of the graph that yields the minimum cost and 
//has the specified divergence values at every node.
//
//To compile: g++ mincost.cc -o mincost
//To run: ./mincost file.dat
//
//Make sure that file.dat is formatted correctly as below
//Problem format: Row 1 indicates |N| and |A|. Rows 2 to |A|+1 correspond
// to arcs, indicating the start node, end node, cost, upper capacity of
// each arc. (lower capacity is assumed to be zero.) 
//The last |N| rows correspond to nodes, indicating the supply at each node.
//
//Four example data sets are included and are posted
//by Paul Tseng on his website
//
//Details given about these data sets are:
//mcnf1.dat: |N|=4, |A|=5, optimal cost = 14
//mcnf2.dat: |N|=20, |A|=80, optimal cost = 30464
//mcnf3.dat: |N|=49, |A|=520, optimal cost = 1.73566449E+11 (corrected) (caution: this network has some parallel arcs)
//mcnf4.dat: |N|=10000, |A|=60000, optimal cost = 21514702 (maybe)
//
//
//Also note that in C++ arrays start at index 0 
//(but the data format lists node labels starting from 1),
//so adjustments had to be made accordingly as to how the data was processed
//In particular, we will subtract 1 from each node label 
//(i.e. node 1 is instead referred to as node 0).
//
//--------------------------------------------------------//
//Algorithm Description:
//
//Complimentary Epsilon Slackness Condition:
//
//A given flow x and potential u satisfy epsilon-complimentary slackness if
//for all arcs j~(i,i'), d[j] + u[i] - u[i'] is
//1. in the interval [-epsilon, epsilon] if 0 < x[j] < c[j]
//2. >= -epsilon if 0 = x[j]
//3. <= epsilon if x[j] = c[j]
//
//Algorithm steps:
//
//Initialization: 
//1. Set u[i] = 0 for all nodes i.
//2. For each arc j~(i,i'), compute r[j] = d[j] + u[i] - u[i']
//3. Initialize flow x by calling initflow(x,r,c,A)
//4. Set e = epsilon = maximum degree of any node so that x and u
//   satisfy epsilon-complimentary slackness condition.
//5. Compute node surplus s[i]= (net flow into node i) - b[i]
//
//Main algorithm:
//1. If s[i] = 0 for all nodes i, then x is feasible, so stop.
//   Else pick any node ibar with s[ibar] > 0.
//   
//   while( s[ibar] != 0 )
//   1. paint each arc j~(i,i'):
//      - white if r[j] \in [-e, -e/2] and x[j] < c[j]
//      - black if r[j] \in [e/2, e] and x[j] > 0
//      - red otherwise
//   2. If there is a black arc jbar out of ibar, then decrease
//        x[jbar] by min{x[jbar], s[ibar]}.
//      else if there is a white arc jbar into ibar, then increase
//        x[jbar] by min{c[jbar] - x[jbar], s[ibar]}.
//   3. Else increase u[ibar] 
//      by as much as possible while maintaining e-comlimentary slackness
//      In particular, set u[ibar] = u[ibar] + alpha where
//      alpha = min{ (-r[j] + e such that j~(ibar, i') and x[j] > 0), 
//                   ( r[j] + e such that j~(i',ibar) and x[j] < c[j]) }
//   end if
//   end while
//   go to 1.
//
//------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <deque>
#include <sys/time.h>

using namespace std;

//given the vector v (for our uses it is a vector of node labels),
//minpot1 finds the minimum value of -r when restricted to the vector v
//of indices (node labels).
double minpot1( vector<int> &v, vector<double> &r ){
  double m;
  m = -r[v.back()];
  v.pop_back();

  while(v.empty() != true){
    m = min(m, -r[v.back()]);
    v.pop_back();
  }
  return m;
}

//same as minpot1, except that minpot2 finds the minimum value component of r 
//(not -r) when restricted to vector v of indices (node labels)
double minpot2( vector<int> &v, vector<double> &r ){
  double m;
  
  m = r[v.back()];
  v.pop_back();

  while(v.empty() != true){
    m = min(m, r[v.back()]);
    v.pop_back();
  }
  return m;
}

//initializes the flow x given a potential u. r is related to u
//by the equation r_j = d_j + u_i - u_i' where j = (i, i') is an arc
void initflow( vector<int> &x, vector<double> const &r, vector<int> const &c, int A){
  int i;
  for(i=0; i < A; i++){
    if( r[i] >= 0){
      x[i] = 0;
    }else{
      x[i] = c[i];
    }
  }
}

//computes the cost of a given flow
double cost( vector<int> const &x, vector<int> const &d, int A){
  int i;
  double totalcost = 0;
  for(i = 0; i < A; i++){
    totalcost += double(x[i])*double(d[i]);
  }
  return totalcost;
}

//checks whether a given flow x is feasible with respect to capacities
//and whether s = 0 (s is a vector indexed by nodes giving the surplus of supply at each node)
bool checkfeas( vector<int> const &x, vector<int> const &c, vector<int> const &s, int A, int N){
  int i;
  for(i=0; i < A; i++){
    if( (x[i] > c[i]) || (x[i]<0) ){ return false; }
  }
  for(i=0; i < N; i++){
    if(s[i] != 0){ return false; }
  }
  return true;
}

void printflow( vector<int> &x, int A){
  int i;
  for(i=0; i < A; i++){
    cout<<"x["<<i<<"] = "<<x[i]<<endl;
  }
}

double diff_sec(timeval t1, timeval t2)
{
    return (t2.tv_sec - t1.tv_sec +  (t2.tv_usec - t1.tv_usec)/1000000.0 );
}


void readdata(const char *filename, int &N, int &A, vector<int> &startnode, vector<int> &endnode, vector<int> &d, vector<int> &c, vector<int> &b, vector<int> &degree){
  ifstream inputdata;
  int databuffer;  //variable storing the last read integer from the file.dat
  int i;
  //Reading in the number of nodes and arcs from the first line of file.dat
  inputdata.open(filename);
  if(!inputdata){
  	cerr << "Can't open input file " << filename << endl;
  	exit(1);
 	}
  inputdata >> N;
  inputdata >> A;
  
  //resize vectors accordingly
  startnode.resize(A);
  endnode.resize(A);
  d.resize(A);
  c.resize(A);
  b.resize(N);
  degree.resize(N);

  //properly import the data into the arrays startnode[], endnode[], c[], d[]
  for(i = 0; i < A; i++){
  	inputdata >> databuffer;
  	startnode[i] = databuffer - 1;
  	degree[databuffer-1] += 1;
  	
  	inputdata >> databuffer;
  	endnode[i] = databuffer - 1;
  	degree[databuffer-1] += 1;
  	
  	inputdata >> d[i];
  	inputdata >> c[i];
  }
  
  //import supply data into b[]
  for(i=0; i<N; i++){
  	inputdata >> b[i];
  }
  inputdata.close();
}


int main( int argc, char *argv[] ) {

  int i; int j;
  int m; int n;
  int N; //number of nodes
  int A; //number of arcs
  const char *filename = argv[1];
  timeval t1;
  timeval t2;
//  time_t t_i; //start time of program
//  time_t t_f; //end time of program
  double dif; //run time of program

  vector<int> startnode;
  vector<int> endnode;
  vector<int> d; //unit cost on each arc
  double d_max; //maximum of absolute values elements in d
  vector<int> c; //upper capacity on each arc
  vector<int> b; //supply at each node
  vector<int> degree; //degree of each node, initially 0  
  int maxdeg;  //max degree of any given node in network
  int maxdegnode; //node index with maximal degree
  
	gettimeofday(&t1, NULL); //records start time of program

  //load data
  readdata(filename, N, A, startnode, endnode, d, c, b, degree);
  
  //Compute the maximum degree and argmax degree over all nodes
  maxdeg=0;
  for(j=0; j<N; j++){ 
    if(maxdeg < degree[j]){ 
      maxdegnode = j;
      maxdeg = degree[j];
    }
  }

  vector< vector<int> > arcout(N, vector<int>(maxdeg));  //Records the indices from the array startnode[] corresponding to arcs out of each node
  vector< vector<int> > arcin(N, vector<int>(maxdeg));  //Records the indices from the array endnode[] corresponding to arcs into each node

  vector<int> numarcout(N,0);
  vector<int> numarcin(N,0);  //initialize so that there are zero arcs coming in or out of each node

  //build arcout and arcin
  for(i=0; i < A; i++){
    m = startnode[i];
    n = endnode[i];
    arcin[n][numarcin[n]] = i;
    numarcin[n] += 1;
    arcout[m][numarcout[m]] = i;
    numarcout[m] += 1;
  }

  //compute the maximum of |d_{j}|
  d_max = fabs(double(d[0]));
  for(i=0; i < A; i++){
    if(fabs(double(d[i]))> d_max){ 
      d_max = fabs(double(d[i])); 
    }
  }

  double epsilon = d_max;  //initial choice of epsilon (large enough to zero node potential satisfies epsilon complimentary slackness
  vector<double> u(N, 0); //node potential, initially zero
  vector<int> x(A, 0); //current flow
  vector<double> r(A); //reduced cost d_j + u_i - u_i' for each arc
  vector<int> s(N); //excess supply at each node, s = Ex - b
  deque<int> nodeq;
  int ibar;

  //Initialize reduced cost r
  for(i=0; i < A; i++){
    r[i] = double(d[i]) + u[startnode[i]] - u[endnode[i]];
  }

  //---------------------------------------------------------------------//
  //Main part of algorithm

  bool proceed;
  int k; //k is always used as a temporary storage variable to simplify the code a bit
  double alpha;
  double alpha1;
  double alpha2;
  int beta;
  vector<int> temp;
  int arc; //temporary variable to store index of a given arc

  while(epsilon >= 1/( double(N) )){
    
    //Change flow so that it satisfies epsilon complementary slackness 
    //with respect to new epsilon.

    initflow(x,r,c,A);

    //Compute s and put nodes with positive s in queue
    for(i=0; i < N; i++){
      s[i] = -b[i];
      for(j=0; j < numarcout[i]; j++){
        k = arcout[i][j];
        s[i] += x[k];
      }
      for(j=0; j < numarcin[i]; j++){
        k = arcin[i][j];
        s[i] -= x[k];
      } 
      if(s[i]>0){
        nodeq.push_back(i);
      }
    }

    while(nodeq.empty() == false){
      ibar = nodeq.back();
      nodeq.pop_back();
      while(s[ibar] > 0){
      	proceed = true;
				temp.clear();

				for(j=0; j<numarcout[ibar] && proceed == true; j++){
					k = arcout[ibar][j];
					if( (r[k] <= epsilon) && (r[k] >= (epsilon/2)) && (x[k] > 0) ){ 
		    		arc = k; 
		    		proceed = false;
      	  }
        }

				if(proceed == false ){
					//if proceed == false, then there are black arcs out of ibar
					k = arc;
					beta = min( x[k], s[ibar] );
					x[k] -= beta;
					s[ibar] -= beta;
					s[endnode[k]] += beta;

					if( (s[endnode[k]] > 0) && (s[endnode[k]] - beta <= 0) ){ 
						nodeq.push_front(endnode[k]);
					}
				}else{
					//This portion executes if there are no black arcs out of ibar
					for( j=0; j<numarcin[ibar] && proceed==true; j++ ){
						k = arcin[ibar][j];
						if( (r[k] >= (-epsilon)) && (r[k] <= (-epsilon/2)) && (x[k] < c[k]) ){
							arc = k;
							proceed = false;
							//Checks for white arcs into ibar. 
							//Loop stops after first arc is found
						}
					}
				
					if(proceed == false ){
						k = arc;
						beta = min( c[k] - x[k], s[ibar]);
						x[k] += beta;
						s[ibar] -= beta;
						s[startnode[k]] += beta;

						if( (s[startnode[k]] > 0) && (s[startnode[k]] - beta) <= 0 ){
							nodeq.push_front(startnode[k]);
						}
					}
				}
				//Increase potential u[ibar] if all arcs are red
				if( proceed==true ){  
					//the vector temp will first store the labels of all arcs out of
					//ibar with positive flow
					for(j=0; j < numarcout[ibar]; j++){
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

						for(j=0; j < numarcin[ibar]; j++){
							//Now we use temp to store the labels of all arcs going
							//into ibar with flow strictly less than the upper capacity.
							k = arcin[ibar][j];
							if( x[k] < c[k] ){ temp.push_back(k); }
						}
						if(temp.empty() == false){
							alpha2 = minpot2( temp, r );
							temp.clear();
							alpha = min(alpha1, alpha2) + epsilon;
							//minpot2(temp, r) minimizes r[j] over all arcs j in temp,
							//which are arcs into ibar with flow strictly less than
							//the upper capacity.
							//We then take alpha to be the minimum of alpha1 and alpha2 + e.
						}else{
							alpha = alpha1 + epsilon;
							//in the case that there are no arcs into ibar with
							//flow strictly less than the upper capacity.
						}
					}else{
						for(j=0; j < numarcin[ibar]; j++){
						k = arcin[ibar][j];
						if( x[k] < c[k] ){ temp.push_back(k); }
						}
						alpha = minpot2( temp, r) + epsilon;
						temp.clear();
					}

				  //update u
				  u[ibar] = u[ibar] + alpha;

				  //update r
				  for(j=0; j < numarcout[ibar]; j++){
				    k = arcout[ibar][j];
				    r[k] += alpha;
				  }
				  for(j=0; j < numarcin[ibar]; j++){
				    k = arcin[ibar][j];
				    r[k] -= alpha;
				  }
    		}     
  		}
		}
		//scale epsilon down by 1/2
		epsilon = epsilon/2;
	}

  cout.setf(ios::floatfield);
  cout.precision(15); //sets number of digits to print out
  
  gettimeofday(&t2, NULL); //records end time of program
  dif = diff_sec(t1, t2);
  cout<<"Run time: "<<dif<<" seconds"<<endl;
  
  //return cost associated with obtained flow x
  if( checkfeas(x,c,s,A,N) ){ 
    cout<<"Current flow is feasible"<<endl;
    cout<<"Minimum cost is: "<<cost(x,d,A)<<endl;
  }else{
    cout<<"Current flow not feasible"<<endl;
  }
}
        
      
       
