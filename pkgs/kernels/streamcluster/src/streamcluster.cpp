/*
 * Copyright (C) 2008 Princeton University
 * All rights reserved.
 * Authors: Jia Deng, Gilberto Contreras
 *
 * streamcluster - Online clustering algorithm
 *
 */
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <sys/resource.h>
#include <climits>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <vector>

#ifdef ENABLE_THREADS
#include <pthread.h>
#include "parsec_barrier.hpp"
#endif

#ifdef TBB_VERSION
#define TBB_STEALER (tbb::task_scheduler_init::occ_stealer)
#define NUM_DIVISIONS (nproc)
#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/cache_aligned_allocator.h"
using namespace tbb;
#endif

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

using namespace std;

#define MAXNAMESIZE 1024 // max filename length
#define SEED 1
/* increase this to reduce probability of random error */
/* increasing it also ups running time of "speedy" part of the code */
/* SP = 1 seems to be fine */
#define SP 1 // number of repetitions of speedy must be >=1

/* higher ITER --> more likely to get correct # of centers */
/* higher ITER also scales the running time almost linearly */
#define ITER 3 // iterate ITER* k log k times; ITER >= 1

#define CACHE_LINE 32 // cache line in byte

/* this structure represents a point */
/* these will be passed around to avoid copying coordinates */
struct Point{
  float weight;
  vector<float>::iterator coord;
  long assign;  /* number of point where this one is assigned */
  float cost;  /* cost of that assignment, weight*distance */
};

/* this is the array of points */
struct Points{
  long num; /* number of points; may not be N if this is a sample */
  int dim;  /* dimensionality */
  vector<Point> p; /* the array itself */
};

static vector<bool> switch_membership; //whether to switch membership in pgain
static vector<bool> is_center; //whether a point is a center
static vector<int> center_table; //index table of centers

static int nproc; //# of threads


#ifdef TBB_VERSION
tbb::cache_aligned_allocator<float> memoryFloat;
tbb::cache_aligned_allocator<Point> memoryPoint;
tbb::cache_aligned_allocator<long> memoryLong;
tbb::cache_aligned_allocator<int> memoryInt;
tbb::cache_aligned_allocator<bool> memoryBool;
#endif


float dist(Point p1, Point p2, long dim);


#ifdef TBB_VERSION
struct HizReduction {
private:
  double hiz;
public:
  Points *points;
  HizReduction(Points *points_): hiz(0),points(points_){}
  HizReduction( HizReduction &d, tbb::split){hiz=0; points = d.points;}

  void operator()(const tbb::blocked_range<int>& range) {
    double myhiz = 0;
    long ptDimension = points->dim;
    int begin = range.begin();
    int end = range.end();
    
    for(int kk=begin; kk!=end; kk++) {
      myhiz += dist(points->p[kk], points->p[0],
			 ptDimension)*points->p[kk].weight;
    }
    hiz += myhiz;
  }

  void join(HizReduction &d){hiz += d.getHiz(); /*fprintf(stderr,"reducing: %lf\n",hiz);*/}
  double getHiz(){return hiz;}

};


struct CenterCreate {
  Points *points;
  CenterCreate(Points *p): points(p){}
  void operator()(const tbb::blocked_range<int>&range) const {
    int begin = range.begin();
    int end = range.end();
    
     for( int k = begin; k!=end; k++ )    {
       float distance = dist(points->p[k],points->p[0],points->dim);
       points->p[k].cost = distance * points->p[k].weight;
       points->p[k].assign=0;
     } 
  }

};



struct CenterOpen {
private:
  double total_cost;
public:
  Points *points;
  int i;
  int type; /*type=0: compute. type=1: reduction */
  CenterOpen(Points *p):points(p),total_cost(0),type(0){}
  CenterOpen(CenterOpen &rhs, tbb::split) 
  {
    total_cost = 0; 
    points = rhs.points;
    i = rhs.i;
    type = rhs.type;
  }

  void operator()(const tbb::blocked_range<int> &range) {
    int begin = range.begin();
    int end = range.end();

    if(type) {
      double local_total = 0.0;
      for(int k = begin; k!=end; k++ )  
	local_total+=points->p[k].cost;
      total_cost += local_total;
    }
    else {
      for(int k = begin; k!=end; k++ )  {
	float distance = dist(points->p[i],points->p[k],points->dim);
	if( i && distance*points->p[k].weight < points->p[k].cost )  {
	  points->p[k].cost = distance * points->p[k].weight;
	  points->p[k].assign=i;
	}
      }
    }
    
  }


  void join(CenterOpen &lhs){total_cost+=lhs.getTotalCost();}
  double getTotalCost(){return total_cost;}

};



class CenterTableCount: public tbb::task{
private:
  Points *points;
  double *work_mem;
  int stride;
  int pid;
public:
  CenterTableCount(int id, int s, Points *p, double *mem):
    pid(id), stride(s), points(p),work_mem(mem){}

  task *execute() {
    int count = 0;
    long bsize = points->num/((NUM_DIVISIONS));
    long k1 = bsize * pid;
    long k2 = k1 + bsize;

    if( pid == (NUM_DIVISIONS)-1 ) 
      k2 = points->num;

    /* fprintf(stderr,"\t[CenterTableCount]: pid=%d stride=%d from %d to %d\n",
       pid, stride, k1, k2); */

    for( int i = k1; i < k2; i++ ) {
      if( is_center[i] ) {
	center_table[i] = count++;
      }
    }

    work_mem[pid*stride] = count;
    //fprintf(stderr,"PID %d done!\n",pid);
    return NULL;
  }

};


class CenterTableCountTask: public tbb::task {
  int is_continuation;
  Points *points;
  double *work_mem;
  int stride;
public:
  CenterTableCountTask(int s, Points *p, double *mem):
    stride(s), points(p), work_mem(mem), is_continuation(0){} 

  task *execute() {
    tbb::task_list list;
    int p;
    
    if(!is_continuation) {
      recycle_as_continuation();
      set_ref_count(NUM_DIVISIONS);

      for(p = 1; p < (NUM_DIVISIONS); p++ ) 
	  list.push_back( *new( allocate_child() ) CenterTableCount(p, stride, points, work_mem));
      CenterTableCount &me = *new( allocate_child() ) CenterTableCount(0, stride, points, work_mem);
      spawn(list);
      is_continuation = 1;
      
      return &me;

    }else {
      /* continuation part */
      int accum = 0;
      for( int p = 0; p < (NUM_DIVISIONS); p++ ) {
	int tmp = (int)work_mem[p*stride];
	work_mem[p*stride] = accum;
	accum += tmp;
      }
      //fprintf(stderr,"Accum = %d\n",accum);
      return NULL;
    }
  }
};


class FixCenter: public tbb::task {
  Points *points;
  double *work_mem;
  int pid;
  int stride;
public:
  FixCenter(int id, int s, Points *p, double *mem):
    pid(id),stride(s),points(p),work_mem(mem){}
  task *execute(){
#ifdef SERIAL_FIXCENTER
    long k1 = 0;
    long k2 = points->num;
#else    
    long bsize = points->num/((NUM_DIVISIONS));
    long k1 = bsize * pid;
    long k2 = k1 + bsize;
    if( pid == (NUM_DIVISIONS)-1 ) k2 = points->num;
#endif
    /*fprintf(stderr,"\t[FixCenter]: pid=%d stride=%d from %d to %d is_center=0x%08x\n",
      pid, stride, k1, k2,(int)is_center);  */
    
    for( int i = k1; i < k2; i++ ) {
      if( is_center[i] ) {
	center_table[i] += (int)work_mem[pid*stride];
	//fprintf(stderr,"\tcenter_table[%d] = %d\n",i,center_table[i]);
      }

    }
      //fprintf(stderr,"PID %d done!\n",pid);
    return NULL;

  }
};

class FixCenterTask: public tbb::task {
  bool is_continuation;
  Points *points;
  double *work_mem;
  int stride;
public:
  FixCenterTask(int s, Points *p, double *mem):
    stride(s), points(p), work_mem(mem), is_continuation(false){} 

  task *execute() {
    tbb::task_list list;
    int p;
    if(!is_continuation) {
      recycle_as_continuation();
      set_ref_count(NUM_DIVISIONS);
      for(p = 1; p < (NUM_DIVISIONS); p++ ) 
	  list.push_back( *new( allocate_child() ) FixCenter(p, stride, points, work_mem));
      spawn(list);
      FixCenter &me = *new (allocate_child()) FixCenter(0, stride, points, work_mem);
      is_continuation = true;
      return &me;
    }else {
      /* coninuation */
      return NULL;
    }
  }
};


class LowerCost: public tbb::task {
  Points *points;
  double *work_mem;
  long x;
  int K;
  int pid;
  int stride;
public:
  LowerCost(int id, int s, Points *p, long x_, double *mem, int k): 
    pid(id), stride(s), points(p), work_mem(mem), K(k), x(x_){}
  task *execute() {

    //my *lower* fields
    double* lower = &work_mem[pid*stride];
    double local_cost_of_opening_x = 0;
    long bsize = points->num/((NUM_DIVISIONS)); //points->num/1;//((NUM_DIVISIONS));
    long k1 = bsize * pid;
    long k2 = k1 + bsize;
    int i;

    if( pid == (NUM_DIVISIONS)-1 ) 
      k2 = points->num;


    /*fprintf(stderr,"\t[LowerCost]: pid=%d stride=%d from %d to %d\n",
      pid, stride, k1, k2);  */
    
    double *cost_of_opening_x = &work_mem[pid*stride + K+1];

    for ( i = k1; i < k2; i++ ) {
      float x_cost = dist(points->p[i], points->p[x], points->dim) 
	* points->p[i].weight;
      float current_cost = points->p[i].cost;

      //fprintf(stderr,"\t (x_cost=%lf < current_cost=%lf)\n",x_cost, current_cost);
      if ( x_cost < current_cost ) {

	// point i would save cost just by switching to x
	// (note that i cannot be a median, 
	// or else dist(p[i], p[x]) would be 0)
	
	switch_membership[i] = 1;
	local_cost_of_opening_x += x_cost - current_cost;
	
      } else {
	
	// cost of assigning i to x is at least current assignment cost of i
	
	// consider the savings that i's **current** median would realize
	// if we reassigned that median and all its members to x;
	// note we've already accounted for the fact that the median
	// would save z by closing; now we have to subtract from the savings
	// the extra cost of reassigning that median and its members 
	int assign = points->p[i].assign;
	lower[center_table[assign]] += current_cost - x_cost;
	//fprintf(stderr,"Lower[%d]=%lf\n",center_table[assign], lower[center_table[assign]]);
      }
    }
    
    *cost_of_opening_x = local_cost_of_opening_x;
    return NULL;
  }
  
  
};
  
class LowerCostTask: public tbb::task {
  bool is_continuation;
  Points *points;
  double *work_mem;
  int K;
  long x;
  int stride;
public:
  LowerCostTask(int s, Points *p, long x_, double *mem, int k): 
    stride(s), points(p), work_mem(mem), K(k), x(x_), is_continuation(false){}

  task *execute() {
    tbb::task_list list;
    int p;
    if(!is_continuation) {
      recycle_as_continuation();
      set_ref_count(NUM_DIVISIONS);
      for(p = 1; p < (NUM_DIVISIONS); p++ ) 
	  list.push_back( *new( allocate_child() )  LowerCost(p, stride, points, x, work_mem, K));
      spawn(list);
      LowerCost &me = *new (allocate_child())  LowerCost(0, stride, points, x, work_mem, K);
      is_continuation = true;
      return &me;
    }else {
      /* continuation */
      return NULL;
    }
  }
};




class CenterClose: public tbb::task {
  Points *points;
  double *work_mem;
  double *number_of_centers_to_close;
  double z;
  int pid, stride;
  int K;

public:
  CenterClose(int id, int s, Points *p, double *mem, int k, double z_): 
    pid(id),stride(s),points(p),work_mem(mem),K(k), z(z_){}

  task *execute() {
    double* gl_lower = &work_mem[(NUM_DIVISIONS)*stride];
    double *cost_of_opening_x;
    int local_number_of_centers_to_close = 0;
    long bsize = points->num/((NUM_DIVISIONS)); //
    long k1 = bsize * pid;
    long k2 = k1 + bsize;

    if( pid == (NUM_DIVISIONS)-1 ) 
      k2 = points->num;

    /*fprintf(stderr,"\t[CenterClose]: pid=%d stride=%d from %d to %d\n",
      pid, stride, k1, k2); */

    number_of_centers_to_close = &work_mem[pid*stride + K];
    cost_of_opening_x = &work_mem[pid*stride + K+1];
    
      for ( int i = k1; i < k2; i++ ) {
	if( is_center[i] ) {
	  double low = z;
	  //aggregate from all threads
	  for( int p = 0; p < (NUM_DIVISIONS); p++ ) {
	    low += work_mem[center_table[i]+p*stride];
	  }
	  gl_lower[center_table[i]] = low;
	  if ( low > 0 ) {
	    // i is a median, and
	    // if we were to open x (which we still may not) we'd close i
	    
	    // note, we'll ignore the following quantity unless we do open x
	    ++local_number_of_centers_to_close;  
	    *cost_of_opening_x -= low;
	  }
	}
      }
      *number_of_centers_to_close = (double)local_number_of_centers_to_close;
      return NULL;
  }

};


class CenterCloseTask: public tbb::task {
  bool is_continuation;
  Points *points;
  double *work_mem;
  int stride;
  double z;
  int K;
public:
  CenterCloseTask(int s, Points *p, double *mem, int k, double z_): 
    stride(s),points(p),work_mem(mem),K(k), z(z_), is_continuation(false){}

  task *execute() {
    tbb::task_list list;
    int p;
    if(!is_continuation) {
      recycle_as_continuation();
      set_ref_count(NUM_DIVISIONS);
      for(p = 1; p < (NUM_DIVISIONS); p++ ) 
	list.push_back( *new( allocate_child() )  CenterClose(p, stride, points, work_mem, K, z));
      spawn(list);
      CenterClose &me = *new (allocate_child())  CenterClose(0, stride, points, work_mem, K, z);
      is_continuation = true;
      return &me;
    }else {
      /* coninuation */


      return NULL;
    }
  }
};



class SaveMoney: public tbb::task{
  Points *points;
  double *work_mem;
  long x;
  int pid, stride;
public:
  SaveMoney(int id, int s, Points *p, long x_, double *mem): 
    pid(id), stride(s), points(p), x(x_), work_mem(mem){}
  task *execute() {
    double* gl_lower = &work_mem[(NUM_DIVISIONS)*stride];
    long bsize = points->num/((NUM_DIVISIONS));//points->num/1;//((NUM_DIVISIONS));
    long k1 = bsize * pid;
    long k2 = k1 + bsize;
    int i;
    
    if( pid == (NUM_DIVISIONS)-1 ) 
      k2 = points->num;

    /*fprintf(stderr,"\t[SaveMoney]: pid=%d stride=%d from %d to %d\n",
      pid, stride, k1, k2);   */
    

    //  we'd save money by opening x; we'll do it
    for ( int i = k1; i < k2; i++ ) {
      bool close_center = gl_lower[center_table[points->p[i].assign]] > 0 ;
      if ( switch_membership[i] || close_center ) {
	// Either i's median (which may be i itself) is closing,
	// or i is closer to x than to its current median
	points->p[i].cost = points->p[i].weight *
	  dist(points->p[i], points->p[x], points->dim);
	points->p[i].assign = x;
	//fprintf(stderr,"\t[SaveMoney] %d: cost %lf, x=%d\n",i,points->p[i].cost, x);
      }
    }
    for( int i = k1; i < k2; i++ ) {
      if( is_center[i] && gl_lower[center_table[i]] > 0 ) {
	is_center[i] = false;
      }
    }
    if( x >= k1 && x < k2 ) {
      //fprintf(stderr,"\t-->is_center[%d]=true!\n",x);
      is_center[x] = true;
    }


    return NULL;
  }
};


class SaveMoneyTask: public tbb::task {
  bool is_continuation;
  Points *points;
  long x;
  double* work_mem;
  int stride;

public:
  SaveMoneyTask(int s, Points *p, long x_, double *mem): 
    stride(s), points(p), x(x_), work_mem(mem) ,is_continuation(false){}


  task *execute() {
    tbb::task_list list;
    int p;
    if(!is_continuation) {
      recycle_as_continuation();
      set_ref_count(NUM_DIVISIONS);
      for(p = 1; p < (NUM_DIVISIONS); p++ ) 
	list.push_back( *new( allocate_child() )  SaveMoney(p, stride, points, x, work_mem));
      spawn(list);
      SaveMoney &me = *new (allocate_child())  SaveMoney(0, stride, points, x, work_mem);
      is_continuation = true;
      return &me;
    }else {
      /* coninuation */


      return NULL;
    }
  }
};

#endif //TBB_VERSION
/********************************************/


/* shuffle points into random order */
void shuffle(Points *points)
{
  for (long i=0;i<points->num-1;i++) {
    long j=(lrand48()%(points->num - i)) + i;
    Point temp = points->p[i];
    points->p[i] = points->p[j];
    points->p[j] = temp;
  }
}

/* shuffle an array of integers */
void intshuffle(vector<long> intarray, long length)
{
  for (long i=0;i<length;i++) {
    long j=(lrand48()%(length - i))+i;
    long temp = intarray[i];
    intarray[i]=intarray[j];
    intarray[j]=temp;
  }
}

/* compute Euclidean distance squared between two points */
float dist(Point p1, Point p2, long dim)
{
  float result=0.0;
  for (int i=0;i<dim;i++){
    result += (p1.coord[i] - p2.coord[i])*(p1.coord[i] - p2.coord[i]);
  }
  return(result);
}

#ifdef TBB_VERSION
/* run speedy on the points, return total cost of solution */
float pspeedy(Points *points, float z, long *kcenter)
{
  static double totalcost;
  static bool open = false;
  static double* costs; //cost for each thread. 
  static int i;


  /* create center at first point, send it to itself */
  {
    int grain_size = points->num / ((NUM_DIVISIONS));
    CenterCreate c(points);
    tbb::parallel_for(tbb::blocked_range<int>(0,points->num, grain_size),c);
  }
    
  *kcenter = 1;


  {
    int grain_size = points->num / ((NUM_DIVISIONS));
    double acc_cost = 0.0;
    CenterOpen c(points);
    for(i = 1; i < points->num; i++ )  {
      bool to_open = ((float)lrand48()/(float)INT_MAX)<(points->p[i].cost/z);
      if( to_open )  {
	(*kcenter)++;
	c.i = i;
	//fprintf(stderr,"** New center for i=%d\n",i);
	tbb::parallel_reduce(tbb::blocked_range<int>(0,points->num,grain_size),c);
      }
    }

    c.type = 1; /* Once last time for actual reduction */
    tbb::parallel_reduce(tbb::blocked_range<int>(0,points->num,grain_size),c);


    totalcost =z*(*kcenter);
    totalcost += c.getTotalCost();
  }
  return(totalcost);
}

#else //!TBB_VERSION

double pspeedy(Points *points, double z, long *kcenter, int pid, [[maybe_unused]] pthread_barrier_t* barrier)
{
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  //my block
  long bsize = points->num/nproc;
  long k1 = bsize * pid;
  long k2 = k1 + bsize;
  if( pid == nproc-1 ){
    k2 = points->num;
  } 

  static double totalcost;

  static bool open = false;
  static vector <double> costs; //cost for each thread. 
  static int i;

#ifdef ENABLE_THREADS
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
#endif

  /* create center at first point, send it to itself */
  for( long k = k1; k < k2; k++ )    {
    float distance = dist(points->p[k],points->p[0],points->dim);
    points->p[k].cost = distance * points->p[k].weight;
    points->p[k].assign=0;
  }

  if( pid==0 )   {
    *kcenter = 1;
    costs.resize(nproc);
  }

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
    
  if( pid != 0 ) { // we are not the master threads. we wait until a center is opened.
    while(true) {
#ifdef ENABLE_THREADS
      pthread_mutex_lock(&mutex);
      while(!open) pthread_cond_wait(&cond,&mutex);
      pthread_mutex_unlock(&mutex);
#endif
      if( i >= points->num ){
        break;
      } 
      for( long k = k1; k < k2; k++ )
	{
	  float distance = dist(points->p[i],points->p[k],points->dim);
	  if( distance*points->p[k].weight < points->p[k].cost )
	    {
	      points->p[k].cost = distance * points->p[k].weight;
	      points->p[k].assign=i;
	    }
	}
#ifdef ENABLE_THREADS
      pthread_barrier_wait(barrier);
      pthread_barrier_wait(barrier);
#endif
    } 
  }
  else  { // I am the master thread. I decide whether to open a center and notify others if so. 
    for(i = 1; i < points->num; i++ )  {
      bool to_open = ((float)lrand48()/(float)INT_MAX)<(points->p[i].cost/z);
      if( to_open )  {
	(*kcenter)++;
#ifdef ENABLE_THREADS
	pthread_mutex_lock(&mutex);
#endif
	open = true;
#ifdef ENABLE_THREADS
	pthread_mutex_unlock(&mutex);
	pthread_cond_broadcast(&cond);
#endif
	for( long k = k1; k < k2; k++ )  {
	  float distance = dist(points->p[i],points->p[k],points->dim);
	  if( distance*points->p[k].weight < points->p[k].cost )  {
	    points->p[k].cost = distance * points->p[k].weight;
	    points->p[k].assign=i;
	  }
	}
#ifdef ENABLE_THREADS
	pthread_barrier_wait(barrier);
#endif
	open = false;
#ifdef ENABLE_THREADS
	pthread_barrier_wait(barrier);
#endif
      }
    }
#ifdef ENABLE_THREADS
    pthread_mutex_lock(&mutex);
#endif
    open = true;
#ifdef ENABLE_THREADS
    pthread_mutex_unlock(&mutex);
    pthread_cond_broadcast(&cond);
#endif
  }
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  open = false;
  double mytotal = 0;
  for( long k = k1; k < k2; k++ )  {
    mytotal += points->p[k].cost;
  }
  costs[pid] = mytotal;
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  // aggregate costs from each thread
  if( pid == 0 )
    {
      totalcost=z*((double)*kcenter);
      for( int i = 0; i < nproc; i++ )
	{
	  totalcost += costs[i];
	} 
    }
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  return(totalcost);
}

#endif // TBB_VERSION


/* For a given point x, find the cost of the following operation:
 * -- open a facility at x if there isn't already one there,
 * -- for points y such that the assignment distance of y exceeds dist(y, x),
 *    make y a member of x,
 * -- for facilities y such that reassigning y and all its members to x 
 *    would save cost, realize this closing and reassignment.
 * 
 * If the cost of this operation is negative (i.e., if this entire operation
 * saves cost), perform this operation and return the amount of cost saved;
 * otherwise, do nothing.
 */

/* numcenters will be updated to reflect the new number of centers */
/* z is the facility cost, x is the number of this point in the array 
   points */


#ifdef TBB_VERSION
double pgain(long x, Points *points, double z, long int *numcenters)
{
  int i;
  int number_of_centers_to_close = 0;

  static double *work_mem;
  static double gl_cost_of_opening_x;
  static int gl_number_of_centers_to_close;

  //each thread takes a block of working_mem.
  int stride = *numcenters+2;

  //make stride a multiple of CACHE_LINE
  int cl = CACHE_LINE/sizeof(double);
  if( stride % cl != 0 ) { 
    stride = cl * ( stride / cl + 1);
  }
  int K = stride -2 ; // K==*numcenters
  
  //my own cost of opening x
  double cost_of_opening_x = 0;

  work_mem = (double*) calloc(stride*((NUM_DIVISIONS)+1),sizeof(double));
  
  gl_cost_of_opening_x = 0;
  gl_number_of_centers_to_close = 0;


  /*For each center, we have a *lower* field that indicates 
    how much we will save by closing the center. 
    Each thread has its own copy of the *lower* fields as an array.
    We first build a table to index the positions of the *lower* fields. 
  */

  /*****  loopA() *****/
  {
    CenterTableCountTask &t = *new ( tbb::task::allocate_root() ) CenterTableCountTask(stride, points, work_mem);
    tbb::task::spawn_root_and_wait(t);
  }

  
  {
    FixCenterTask &t = *new ( tbb::task::allocate_root() ) FixCenterTask(stride, points, work_mem);
    tbb::task::spawn_root_and_wait(t);
  }    

  /***************/

  //now we finish building the table. clear the working memory.
  memset(switch_membership, 0, points->num*sizeof(bool));
  memset(work_mem, 0, (NUM_DIVISIONS+1)*stride*sizeof(double));

  /* loopB */
  {
    LowerCostTask &t = *new ( tbb::task::allocate_root() )  LowerCostTask(stride, points, x, work_mem, K);
    tbb::task::spawn_root_and_wait(t);
  }    

  /* LoopC */
  {
    CenterCloseTask &t = *new ( tbb::task::allocate_root() )  CenterCloseTask(stride, points, work_mem, K, z);
    tbb::task::spawn_root_and_wait(t);
  }    


  gl_cost_of_opening_x = z;
  //aggregate
  for( int p = 0; p < (NUM_DIVISIONS); p++ ) {
    gl_number_of_centers_to_close += (int)work_mem[p*stride + K];
    gl_cost_of_opening_x += work_mem[p*stride+K+1];
  }

  /*fprintf(stderr,"\tgl_number_of_centers_to_close = %d\n",gl_number_of_centers_to_close);
    fprintf(stderr,"\tgl_cost_of_opening_x = %lf\n",gl_cost_of_opening_x); */


  // Now, check whether opening x would save cost; if so, do it, and
  // otherwise do nothing

  if ( gl_cost_of_opening_x < 0 ) {

    /* loopD */
    SaveMoneyTask &t = *new ( tbb::task::allocate_root() )  SaveMoneyTask(stride, points, x, work_mem);
    tbb::task::spawn_root_and_wait(t);


    *numcenters = *numcenters + 1 - gl_number_of_centers_to_close;    
  }
  else {
    gl_cost_of_opening_x = 0;  // the value we'll return
  }

  free(work_mem);

  return -gl_cost_of_opening_x;
}

#else //!TBB_VERSION


double pgain(long x, Points *points, double z, long int *numcenters, int pid, [[maybe_unused]] pthread_barrier_t* barrier)
{
  //  printf("pgain pthread %d begin\n",pid);
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  //my block
  long bsize = points->num/nproc;
  long k1 = bsize * pid;
  long k2 = k1 + bsize;
  if( pid == nproc-1 ){
    k2 = points->num;
  } 

  
  int number_of_centers_to_close = 0;

  static vector<double> work_mem;
  static double gl_cost_of_opening_x;
  static int gl_number_of_centers_to_close;

  //each thread takes a block of working_mem.
  long stride = *numcenters+2;
  //make stride a multiple of CACHE_LINE
  int cl = CACHE_LINE/sizeof(double);
  if( stride % cl != 0 ) { 
    stride = cl * ( stride / cl + 1);
  }
  long K = stride -2 ; // K==*numcenters
  
  //my own cost of opening x
  double cost_of_opening_x = 0;

  if( pid==0 )    { 
    work_mem.resize(stride*(nproc+1));
    gl_cost_of_opening_x = 0;
    gl_number_of_centers_to_close = 0;
  }

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  /*For each center, we have a *lower* field that indicates 
    how much we will save by closing the center. 
    Each thread has its own copy of the *lower* fields as an array.
    We first build a table to index the positions of the *lower* fields. 
  */

  int count = 0;
  for( long i = k1; i < k2; i++ ) {
    if( is_center[i] ) {
      center_table[i] = count++;
    }
  }
  work_mem[pid*stride] = count;

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  if( pid == 0 ) {
    int accum = 0;
    for( int p = 0; p < nproc; p++ ) {
      int tmp = (int)work_mem[p*stride];
      work_mem[p*stride] = accum;
      accum += tmp;
    }
  }

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  for(long i = k1; i < k2; i++ ) {
    if( is_center[i] ) {
      center_table[i] += (int)work_mem[pid*stride];
    }
  }

  //now we finish building the table. clear the working memory.
  fill(switch_membership.begin()+k1, switch_membership.end(), false);
  fill(work_mem.begin()+pid*stride, work_mem.end(), 0);
  if( pid== 0 ){
    fill(work_mem.begin()+nproc*stride, work_mem.end(), 0);
  } 

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  
  //my *lower* fields
  auto lower = work_mem.begin() + pid*stride;
  //global *lower* fields
  auto gl_lower = work_mem.begin() + nproc*stride;

  for (long i = k1; i < k2; i++ ) {
    float x_cost = dist(points->p[i], points->p[x], points->dim) 
      * points->p[i].weight;
    float current_cost = points->p[i].cost;

    if ( x_cost < current_cost ) {

      // point i would save cost just by switching to x
      // (note that i cannot be a median, 
      // or else dist(p[i], p[x]) would be 0)
      
      switch_membership[i] = true;
      cost_of_opening_x += x_cost - current_cost;

    } else {

      // cost of assigning i to x is at least current assignment cost of i

      // consider the savings that i's **current** median would realize
      // if we reassigned that median and all its members to x;
      // note we've already accounted for the fact that the median
      // would save z by closing; now we have to subtract from the savings
      // the extra cost of reassigning that median and its members 
      long assign = points->p[i].assign;
      lower[center_table[assign]] += current_cost - x_cost;
    }
  }

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  // at this time, we can calculate the cost of opening a center
  // at x; if it is negative, we'll go through with opening it

  for ( long i = k1; i < k2; i++ ) {
    if( is_center[i] ) {
      double low = z;
      //aggregate from all threads
      for( int p = 0; p < nproc; p++ ) {
	low += work_mem[center_table[i]+p*stride];
      }
      gl_lower[center_table[i]] = low;
      if ( low > 0 ) {
	// i is a median, and
	// if we were to open x (which we still may not) we'd close i

	// note, we'll ignore the following quantity unless we do open x
	++number_of_centers_to_close;  
	cost_of_opening_x -= low;
      }
    }
  }
  //use the rest of working memory to store the following
  work_mem[pid*stride + K] = number_of_centers_to_close;
  work_mem[pid*stride + K+1] = cost_of_opening_x;

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  //  printf("thread %d cost complete\n",pid); 

  if( pid==0 ) {
    gl_cost_of_opening_x = z;
    //aggregate
    for( int p = 0; p < nproc; p++ ) {
      gl_number_of_centers_to_close += (int)work_mem[p*stride + K];
      gl_cost_of_opening_x += work_mem[p*stride+K+1];
    }
  }
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
  // Now, check whether opening x would save cost; if so, do it, and
  // otherwise do nothing

  if ( gl_cost_of_opening_x < 0 ) {
    //  we'd save money by opening x; we'll do it
    for ( long i = k1; i < k2; i++ ) {
      bool close_center = gl_lower[center_table[points->p[i].assign]] > 0 ;
      if ( switch_membership[i] || close_center ) {
	// Either i's median (which may be i itself) is closing,
	// or i is closer to x than to its current median
	points->p[i].cost = points->p[i].weight *
	  dist(points->p[i], points->p[x], points->dim);
	points->p[i].assign = x;
      }
    }
    for( long i = k1; i < k2; i++ ) {
      if( is_center[i] && gl_lower[center_table[i]] > 0 ) {
	is_center[i] = false;
      }
    }
    if( x >= k1 && x < k2 ) {
      is_center[x] = true;
    }

    if( pid==0 ) {
      *numcenters = *numcenters + 1 - gl_number_of_centers_to_close;
    }
  }
  else {
    if( pid==0 ){
      gl_cost_of_opening_x = 0;  // the value we'll return
    }
  }
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  return -gl_cost_of_opening_x;
}

#endif // TBB_VERSION



/* facility location on the points using local search */
/* z is the facility cost, returns the total cost and # of centers */
/* assumes we are seeded with a reasonable solution */
/* cost should represent this solution's cost */
/* halt if there is < e improvement after iter calls to gain */
/* feasible is an array of numfeasible points which may be centers */

#ifdef TBB_VERSION
float pFL(Points *points, int *feasible, int numfeasible,
	  double z, long *k, double cost, long iter, double e)
{

  long i;
  long x;
  double change;
  long numberOfPoints;

  change = cost;
  /* continue until we run iter iterations without improvement */
  /* stop instead if improvement is less than e */
  while (change/cost > 1.0*e) {
    change = 0.0;
    numberOfPoints = points->num;
    /* randomize order in which centers are considered */    
    intshuffle(feasible, numfeasible);

    for (i=0;i<iter;i++) {
      x = i%numfeasible;
      //fprintf(stderr,"Iteration %d z=%lf, change=%lf\n",i,z,change);
      change += pgain(feasible[x], points, z , k);
      //fprintf(stderr,"*** change: %lf, z=%lf\n",change,z);
    }
    cost -= change;
  }

  return(cost);
}


#else //!TBB_VERSION
 double pFL(Points *points, vector<long> feasible, long numfeasible,
	  double z, long *k, double cost, long iter, float e, 
	  int pid, pthread_barrier_t* barrier)
{
#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif
 

  double change = cost;
  /* continue until we run iter iterations without improvement */
  /* stop instead if improvement is less than e */
  while (change/cost > 1.0*e) {
    change = 0.0;
    
    /* randomize order in which centers are considered */

    if( pid == 0 ) {
      intshuffle(feasible, numfeasible);
    }
#ifdef ENABLE_THREADS
    pthread_barrier_wait(barrier);
#endif
    for (long i=0;i<iter;i++) {
      long x = i%numfeasible;
      change += pgain(feasible[x], points, z, k, pid, barrier);
    }
    cost -= change;
#ifdef ENABLE_THREADS
    pthread_barrier_wait(barrier);
#endif
  }
  return(cost);
}

#endif // TBB_VERSION



#ifdef TBB_VERSION
long selectfeasible_fast(Points *points, long **feasible, long kmin)
#else
long selectfeasible_fast(Points *points, vector<long> *feasible, long kmin,[[maybe_unused]] int pid, [[maybe_unused]] pthread_barrier_t* barrier)
#endif
{
  long numfeasible = points->num;
  if (numfeasible > (long)(ITER*kmin*log((double)kmin))){
    numfeasible = (long)(ITER*kmin*log((double)kmin));
  }
  feasible->resize(numfeasible);
  

  /* 
     Calcuate my block. 
     For now this routine does not seem to be the bottleneck, so it is not parallelized. 
     When necessary, this can be parallelized by setting k1 and k2 to 
     proper values and calling this routine from all threads ( it is called only
     by thread 0 for now ). 
     Note that when parallelized, the randomization might not be the same and it might
     not be difficult to measure the parallel speed-up for the whole program. 
   */
  //  long bsize = numfeasible;
  long k1 = 0;
  long k2 = numfeasible;

  /* not many points, all will be feasible */
  if (numfeasible == points->num) {
    for (long i=k1;i<k2;i++){
      (*feasible)[i] = i;
    }
    return numfeasible;
  }
#ifdef TBB_VERSION
  vector<float> accumweight (points->num, 0.0);
#else
  vector<float> accumweight (points->num, 0.0);
#endif

  accumweight[0] = points->p[0].weight;
  
  for( int i = 1; i < points->num; i++ ) {
    accumweight[i] = accumweight[i-1] + points->p[i].weight;
  }
  float totalweight=accumweight[points->num-1];

  for(long i=k1; i<k2; i++ ) {
    float w = ((float)lrand48()/(float)INT_MAX)*totalweight;
    //binary search
    long l=0;
    long r=points->num-1;
    if( accumweight[0] > w )  { 
      (*feasible)[i]=0; 
      continue;
    }
    while( l+1 < r ) {
      long k = (l+r)/2;
      if( accumweight[k] > w ) {
	r = k;
      } 
      else {
	l=k;
      }
    }
    (*feasible)[i]=r;
  }

  return numfeasible;
}



#ifdef TBB_VERSION
/* compute approximate kmedian on the points */
float pkmedian(Points *points, long kmin, long kmax, long* kfinal,
	       int pid, pthread_barrier_t* barrier )
{
  int i;
  double cost;
  double lastcost;
  double hiz, loz, z;

  static long k;
  static int *feasible;
  static int numfeasible;
  static double* hizs;


  //  hizs = (double*)calloc(nproc,sizeof(double));
  hiz = loz = 0.0;
  long numberOfPoints = points->num;
  long ptDimension = points->dim;

  //my block
  long bsize = points->num/nproc;
  long k1 = bsize * pid;
  long k2 = k1 + bsize;
  if( pid == nproc-1 ) k2 = points->num;

  
  //fprintf(stderr,"Starting Kmedian procedure\n");
  //fprintf(stderr,"%i points in %i dimensions\n", numberOfPoints, ptDimension);

  int grain_size = points->num / ((NUM_DIVISIONS));
  if(grain_size==0)
    {
      
      for (long kk=0;kk < points->num; kk++ ) 
	{
	  hiz += dist(points->p[kk], points->p[0],
		      ptDimension)*points->p[kk].weight;
	}
      
    }
  else {
    HizReduction h(points);
    tbb::parallel_reduce(tbb::blocked_range<int>(0,points->num, grain_size), h);
    hiz = h.getHiz();
  }

  loz=0.0; z = (hiz+loz)/2.0;

  /* NEW: Check whether more centers than points! */
  if (points->num <= kmax) {
    /* just return all points as facilities */
      for (long kk=0;kk<points->num;kk++) 
	{
	  points->p[kk].assign = kk;
	  points->p[kk].cost = 0;
	}
    
    cost = 0;
    *kfinal = k;

    return cost;
  }

    shuffle(points);
    cost = pspeedy(points, z, &k);

    i=0;

  /* give speedy SP chances to get at least kmin/2 facilities */
  while ((k < kmin)&&(i<SP)) {
    cost = pspeedy(points, z, &k);
    i++;
  }

  /* if still not enough facilities, assume z is too high */
  while (k < kmin) {
    if (i >= SP) 
      {hiz=z; z=(hiz+loz)/2.0; i=0;}
    
    shuffle(points);
    cost =  pspeedy(points, z, &k);
    i++;
  }

  /* now we begin the binary search for real */
  /* must designate some points as feasible centers */
  /* this creates more consistancy between FL runs */
  /* helps to guarantee correct # of centers at the end */

    numfeasible = selectfeasible_fast(points,&feasible,kmin);
    for( int i = 0; i< points->num; i++ ) {
      //fprintf(stderr,"\t-->is_center[%d]=true!\n",points->p[i].assign);
      is_center[points->p[i].assign]= true;
    }


  while(1) {
    /* first get a rough estimate on the FL solution */
    lastcost = cost;
    cost = pFL(points, feasible, numfeasible,
	       z, &k, cost, (long)(ITER*kmax*log((double)kmax)), 0.1);

    /* if number of centers seems good, try a more accurate FL */
    if (((k <= (1.1)*kmax)&&(k >= (0.9)*kmin))||
	((k <= kmax+2)&&(k >= kmin-2))) {
      
      /* may need to run a little longer here before halting without
	 improvement */
      cost = pFL(points, feasible, numfeasible,
		 z, &k, cost, (long)(ITER*kmax*log((double)kmax)), 0.001);
    }

    if (k > kmax) {
      /* facilities too cheap */
      /* increase facility cost and up the cost accordingly */
      loz = z; z = (hiz+loz)/2.0;
      cost += (z-loz)*k;
    }
    if (k < kmin) {
      /* facilities too expensive */
      /* decrease facility cost and reduce the cost accordingly */
      hiz = z; z = (hiz+loz)/2.0;
      cost += (z-hiz)*k;
    }

    /* if k is good, return the result */
    /* if we're stuck, just give up and return what we have */
    if (((k <= kmax)&&(k >= kmin))||((loz >= (0.999)*hiz)) )
      { 
	break;
      }

  }

  //  fprintf(stderr,"Cleaning up...\n");
  //clean up...
  free(feasible); 
  *kfinal = k;

  return cost;
}


#else //!TBB_VERSION

/* compute approximate kmedian on the points */
double pkmedian(Points *points, long kmin, long kmax, long* kfinal,
	       int pid, pthread_barrier_t* barrier )
{
  static long k;
  static vector<long> feasible;
  static long numfeasible;
  static vector<double> hizs;

  if( pid==0 ){
    hizs.resize(nproc);
  }
  double hiz = 0.0;

  long ptDimension = points->dim;

  //my block
  long bsize = points->num/nproc;
  long k1 = bsize * pid;
  long k2 = k1 + bsize;
  if( pid == nproc-1 ){
    k2 = points->num;
  } 

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  double myhiz = 0;
  for (long kk=k1;kk < k2; kk++ ) {
    myhiz += dist(points->p[kk], points->p[0],
		      ptDimension)*points->p[kk].weight;
  }
  hizs[pid] = myhiz;

#ifdef ENABLE_THREADS  
  pthread_barrier_wait(barrier);
#endif

  for( int i = 0; i < nproc; i++ )   {
    hiz += hizs[i];
  }

  /*Constexpr to avoid magic numbers */
  constexpr double divisor = 2.0;
  constexpr double pfl_e = 0.1;
  constexpr double pfl_e_smaller = 0.001;
  constexpr double fl_max = 1.1;
  constexpr double fl_min = 0.9;
  constexpr double fl_almost_one = 0.999;

  double loz=0.0; 
  double z = (hiz+loz)/divisor;
  /* NEW: Check whether more centers than points! */
  if (points->num <= kmax) {
    /* just return all points as facilities */
    for (long kk=k1;kk<k2;kk++) {
      points->p[kk].assign = kk;
      points->p[kk].cost = 0;
    }
    double cost = 0;
    if( pid== 0 ) {
      *kfinal = k;
    }
    return cost;
  }

  if( pid == 0 ){
    shuffle(points);
  } 
  double cost = pspeedy(points, z, &k, pid, barrier);

  int i=0;
  /* give speedy SP chances to get at least kmin/2 facilities */
  while ((k < kmin)&&(i<SP)) {
    cost = pspeedy(points, z, &k, pid, barrier);
    i++;
  }

  /* if still not enough facilities, assume z is too high */
  while (k < kmin) {
    if (i >= SP) {hiz=z; z=(hiz+loz)/divisor; i=0;}
    if( pid == 0 ){
      shuffle(points);
    } 
    cost = pspeedy(points, z, &k, pid, barrier);
    i++;
  }

  /* now we begin the binary search for real */
  /* must designate some points as feasible centers */
  /* this creates more consistancy between FL runs */
  /* helps to guarantee correct # of centers at the end */
  
  if( pid == 0 )
    {
      numfeasible = selectfeasible_fast(points,&feasible,kmin,pid,barrier);
      for( int i = 0; i< points->num; i++ ) {
	is_center[points->p[i].assign]= true;
      }
    }

#ifdef ENABLE_THREADS
  pthread_barrier_wait(barrier);
#endif

  while(true) {
    /* first get a rough estimate on the FL solution */
    
    cost = pFL(points, feasible, numfeasible,
	       z, &k, cost, (long)(ITER*kmax*log((double)kmax)), pfl_e, pid, barrier);

    /* if number of centers seems good, try a more accurate FL */
    if ((((double)k <= (fl_max)*(double)kmax)&&((double)k >= (fl_min)*(double)kmin))||
	((k <= kmax+2)&&(k >= kmin-2))) {

      /* may need to run a little longer here before halting without
	 improvement */
      cost = pFL(points, feasible, numfeasible,
		 z, &k, cost, (long)(ITER*kmax*log((double)kmax)), pfl_e_smaller, pid, barrier);
    }

    if (k > kmax) {
      /* facilities too cheap */
      /* increase facility cost and up the cost accordingly */
      loz = z; z = (hiz+loz)/divisor;
      cost += (z-loz)*(double)k;
    }
    if (k < kmin) {
      /* facilities too expensive */
      /* decrease facility cost and reduce the cost accordingly */
      hiz = z; z = (hiz+loz)/divisor;
      cost += (z-hiz)*(double)k;
    }

    /* if k is good, return the result */
    /* if we're stuck, just give up and return what we have */
    if (((k <= kmax)&&(k >= kmin))||((loz >= (fl_almost_one)*hiz)) )
      { 
	break;
      }
#ifdef ENABLE_THREADS
    pthread_barrier_wait(barrier);
#endif
  }

  //clean up...
  if( pid==0 ) {
    *kfinal = k;
  }

  return cost;
}

#endif // TBB_VERSION




/* compute the means for the k clusters */
int contcenters(Points *points)
{

  for (long i=0;i<points->num;i++) {
    /* compute relative weight of this point to the cluster */
    if (points->p[i].assign != i) {
      float relweight=points->p[points->p[i].assign].weight + points->p[i].weight;
      relweight = points->p[i].weight/relweight;
      for (long ii=0;ii<points->dim;ii++) {
	points->p[points->p[i].assign].coord[ii]*=(float)1.0-relweight;
	points->p[points->p[i].assign].coord[ii]+=
	  points->p[i].coord[ii]*relweight;
      }
      points->p[points->p[i].assign].weight += points->p[i].weight;
    }
  }
  
  return 0;
}

/* copy centers from points to centers */
void copycenters(Points *points, Points* centers, vector<long> centerIDs, long offset)
{

  vector<bool> is_a_median (points->num, false);

  /* mark the centers */
  for (long i = 0; i < points->num; i++ ) {
    is_a_median[points->p[i].assign] = true;
  }

  long k=centers->num;

  /* count how many  */
  for (long i = 0; i < points->num; i++ ) {
    if ( is_a_median[i] ) {
      copy_n(points->p[i].coord, points->dim, centers->p[k].coord);
      centers->p[k].weight = points->p[i].weight;
      centerIDs[k] = i + offset;
      k++;
    }
  }

  centers->num = k;

}

struct pkmedian_arg_t
{
  Points* points;
  long kmin;
  long kmax;
  long* kfinal;
  int pid;
  pthread_barrier_t* barrier;
};

void* localSearchSub(void* arg_) {

  auto arg= static_cast<pkmedian_arg_t*> (arg_);
  pkmedian(arg->points,arg->kmin,arg->kmax,arg->kfinal,arg->pid,arg->barrier);

  return nullptr;
}

#ifdef TBB_VERSION
void localSearch( Points* points, long kmin, long kmax, long* kfinal ) {
  pkmedian_arg_t arg;
  arg.points = points;
  arg.kmin = kmin;
  arg.kmax = kmax;
  arg.pid = 0;
  arg.kfinal = kfinal;
  localSearchSub(&arg);
}
#else //!TBB_VERSION

void localSearch( Points* points, long kmin, long kmax, long* kfinal ) {
    pthread_barrier_t barrier;
    vector<pthread_t> threads(nproc);
    vector<pkmedian_arg_t> arg (nproc);

#ifdef ENABLE_THREADS
    pthread_barrier_init(&barrier,NULL,nproc);
#endif
    for( int i = 0; i < nproc; i++ ) {
      arg[i].points = points;
      arg[i].kmin = kmin;
      arg[i].kmax = kmax;
      arg[i].pid = i;
      arg[i].kfinal = kfinal;

      arg[i].barrier = &barrier;
#ifdef ENABLE_THREADS
      pthread_create(&threads[i],NULL,localSearchSub,(void*)&arg[i]);
#else
      localSearchSub(&arg[0]);
#endif
    }

#ifdef ENABLE_THREADS
    for ( int i = 0; i < nproc; i++) {
      pthread_join(threads[i],NULL);
    }
#endif

#ifdef ENABLE_THREADS
    pthread_barrier_destroy(&barrier);
#endif
}
#endif // TBB_VERSION


class PStream {
public:
  virtual size_t read( float* dest, int dim, int num ) = 0;
  virtual int ferror() = 0;
  virtual int feof() = 0;
  virtual ~PStream() = default;
  PStream() = default;
  PStream(const PStream& other) = default;
  PStream& operator=(const PStream& other) = default;
  PStream(PStream&& other) = default;
  PStream& operator=(PStream&& other) = default;
};

//synthetic stream
class SimStream : public PStream {
public:
  SimStream(long n_ ) {
    n = n_;
  }
  size_t read( float* dest, int dim, int num ) override{
    size_t count = 0;
    for( int i = 0; i < num && n > 0; i++ ) {
      for( int k = 0; k < dim; k++ ) {
	dest[i*dim + k] = (float)lrand48()/(float)INT_MAX;
      }
      n--;
      count++;
    }
    return count;
  }
  int ferror() override{
    return 0;
  }
  int feof() override{
    return static_cast<int> (n <= 0);
  }
  ~SimStream() = default;
  SimStream(const SimStream& other) = default;
  SimStream& operator=(const SimStream& other) = default;
  SimStream(SimStream&& other) = default;
  SimStream& operator=(SimStream&& other) = default;
  
private:
  long n;
};

class FileStream : public PStream {
public:
  FileStream(const string & filename ) :
    fp{fopen(filename.c_str(), "rb"), fclose}
    {
      if( fp.get() == nullptr ) {
        fmt::print(stderr,"error opening file {}\n.",filename);
        exit(1);
      }
    }
  size_t read( float* dest, int dim, int num ) override{
    return std::fread(dest, sizeof(float)*dim, num, fp.get()); 
  }
  int ferror() override{
    return std::ferror(fp.get());
  }
  int feof() override{
    return std::feof(fp.get());
  }
  ~FileStream() {
    fmt::print(stderr,"closing file stream\n");
  }
  FileStream(const FileStream& other) = default;
  FileStream& operator=(const FileStream& other) = default;
  FileStream(FileStream&& other) = default;
  FileStream& operator=(FileStream&& other) = default;
private:
  unique_ptr<FILE, int(*)(FILE*)> fp;
};

void outcenterIDs( Points* centers, vector<long> centerIDs, const string &outfile ) {
  unique_ptr<FILE, int(*)(FILE*)> fp {fopen(outfile.c_str(), "w"), fclose};
  if( fp.get()==nullptr ) {
    fmt::print(stderr, "error opening {}\n",outfile);
    exit(1);
  }
  vector<int> is_a_median(centers->num, 0);
  for( int i =0 ; i< centers->num; i++ ) {
    is_a_median[centers->p[i].assign] = 1;
  }

  for( int i = 0; i < centers->num; i++ ) {
    if( is_a_median[i] != 0) {
      fmt::print(fp.get(), "{}\n", centerIDs[i]);
      fmt::print(fp.get(), "%{}\n", centers->p[i].weight);
      for( int k = 0; k < centers->dim; k++ ) {
  fmt::print(fp.get(), "{} ", centers->p[i].coord[k]);
      }
      fmt::print(fp.get(),"\n\n");
    }
  } 
}

void streamCluster( PStream* stream, 
		    long kmin, long kmax, int dim,
		    long chunksize, long centersize, const string & outfile )
{

#ifdef TBB_VERSION
  vector<float> block (chunksize*dim, 0.0);
  vector<float> centerBlock(centersize*dim, 0.0);
  vector<long> centerIDs (centersize*dim, 0);
#else
  vector<float> block (chunksize*dim, 0.0);
  vector<float> centerBlock(centersize*dim, 0.0);
  vector<long> centerIDs (centersize*dim, 0);
#endif

  if( block.empty() ) { 
    fmt::print(stderr,"not enough memory for a chunk!\n");
    exit(1);
  }

  Points points;
  points.dim = dim;
  points.num = chunksize;
  points.p.resize(chunksize);

  for( int i = 0; i < chunksize; i++ ) {
    points.p[i].coord = block.begin()+(i*dim);
  }

  Points centers;
  centers.dim = dim;
  centers.p.resize(centersize);
  centers.num = 0;

  for( int i = 0; i< centersize; i++ ) {
    centers.p[i].coord = centerBlock.begin()+(i*dim);
    centers.p[i].weight = 1.0;
  }

  long IDoffset = 0;
  long kfinal = 0;
  while(true) {

    float* aux_block = &block[0];
    size_t numRead  = stream->read(aux_block, dim, (int)chunksize ); 
    fmt::print(stderr,"read {} points\n",numRead);

    if( stream->ferror() != 0 || numRead < (unsigned int)chunksize && stream->feof() == 0 ) {
      fmt::print(stderr, "error reading data!\n");
      exit(1);
    }

    points.num = numRead;
    for( int i = 0; i < points.num; i++ ) {
      points.p[i].weight = 1.0;
    }

#ifdef TBB_VERSION
    switch_membership = (bool*)memoryBool.allocate(points.num*sizeof(bool), NULL);
    is_center = (bool*)calloc(points.num,sizeof(bool));
    center_table = (int*)memoryInt.allocate(points.num*sizeof(int));
#else
    switch_membership.resize(points.num); 
    is_center.resize(points.num);
    fill(is_center.begin(), is_center.end(), false);
    center_table.resize(points.num);
#endif


    //fprintf(stderr,"center_table = 0x%08x\n",(int)center_table);
    //fprintf(stderr,"is_center = 0x%08x\n",(int)is_center);

    localSearch(&points,kmin, kmax,&kfinal); // parallel

    //fprintf(stderr,"finish local search\n");
    contcenters(&points); /* sequential */
    if( kfinal + centers.num > centersize ) {
      //here we don't handle the situation where # of centers gets too large. 
      fmt::print(stderr,"oops! no more space for centers\n");
      exit(1);
    }

    copycenters(&points, &centers, centerIDs, IDoffset); /* sequential */
    IDoffset += numRead;

#ifdef TBB_VERSION
    memoryBool.deallocate(switch_membership, sizeof(bool));
    free(is_center);
    memoryInt.deallocate(center_table, sizeof(int));
#endif

    if( stream->feof() != 0) {
      break;
    }
  }

  //finally cluster all temp centers
#ifdef TBB_VERSION
  switch_membership = (bool*)memoryBool.allocate(centers.num*sizeof(bool));
  is_center.resize(centers.num); 
  center_table.resize(centers.num); 
#else
  switch_membership.resize(centers.num);
  is_center.resize(centers.num);
  fill(is_center.begin(), is_center.end(), false);
  center_table.resize(centers.num);
#endif

  localSearch( &centers, kmin, kmax ,&kfinal ); // parallel
  contcenters(&centers);
  outcenterIDs( &centers, centerIDs, outfile);

// Added to solve memory leaks, from here 
#ifdef TBB_VERSION
    memoryPoint.deallocate(points.p, sizeof(Point));
    memoryFloat.deallocate(block, sizeof(float));
    memoryPoint.deallocate(Centers.p, sizeof(Point));
    memoryFloat.deallocate(centerBlock, sizeof(float));
    memoryLong.deallocate(centerIDs, sizeof(long));
#endif
// Up to here

}

int main(int argc, char **argv)
{

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
        fprintf(stderr,"PARSEC Benchmark Suite Version " __PARSEC_XSTRING(PARSEC_VERSION)"\n");
	fflush(NULL);
#else
        fmt::print(stderr,"PARSEC Benchmark Suite\n");
	fflush(nullptr);
#endif //PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
  __parsec_bench_begin(__parsec_streamcluster);
#endif

  constexpr int min_argc = 10;
  
  vector<string> argv_vec (argv, argv+argc);
  if (argc<min_argc) {
    fmt::print(stderr,"usage: {} k1 k2 d n chunksize clustersize infile outfile nproc\n",
	    argv_vec[0]);
    fmt::print(stderr,"  k1:          Min. number of centers allowed\n");
    fmt::print(stderr,"  k2:          Max. number of centers allowed\n");
    fmt::print(stderr,"  d:           Dimension of each data point\n");
    fmt::print(stderr,"  n:           Number of data points\n");
    fmt::print(stderr,"  chunksize:   Number of data points to handle per step\n");
    fmt::print(stderr,"  clustersize: Maximum number of intermediate centers\n");
    fmt::print(stderr,"  infile:      Input file (if n<=0)\n");
    fmt::print(stderr,"  outfile:     Output file\n");
    fmt::print(stderr,"  nproc:       Number of threads to use\n");
    fmt::print(stderr,"\n");
    fmt::print(stderr, "if n > 0, points will be randomly generated instead of reading from infile.\n");
    exit(1);
  }


  int argv_index = 1;
  long kmin = atoi(argv_vec[argv_index++].c_str());
  long kmax = atoi(argv_vec[argv_index++].c_str());
  int dim = atoi(argv_vec[argv_index++].c_str());
  long n = atoi(argv_vec[argv_index++].c_str());
  long chunksize = atoi(argv_vec[argv_index++].c_str());
  long clustersize = atoi(argv_vec[argv_index++].c_str());
  string infilename = argv_vec[argv_index++];
  string outfilename = argv_vec[argv_index++];
  nproc = atoi(argv_vec[argv_index++].c_str());


#ifdef TBB_VERSION
  fprintf(stderr,"TBB version. Number of divisions: %d\n",NUM_DIVISIONS);
  tbb::task_scheduler_init init(nproc);
#endif


  srand48(SEED);
  unique_ptr<PStream> stream;
  if( n > 0 ) {
    stream = make_unique <SimStream> (n);
  }
  else {
    stream = make_unique <FileStream> (infilename); 
  }


#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_begin();
#endif

  streamCluster(stream.get(), kmin, kmax, dim, chunksize, clustersize, outfilename );

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif


#ifdef ENABLE_PARSEC_HOOKS
  __parsec_bench_end();
#endif
  
  return 0;
}
