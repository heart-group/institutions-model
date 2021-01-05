#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <string.h>
#include <math.h>

//Define an agent - an agent can belong to classes S, U, P, R
//Additionally agent can be in contact list
//Additionally agent can get tested and wait for test results
//Time in status and time since test is used by functions to change status
// mi is mean mobility
// ph is the class ratio of high to low mobility
// mh is mobility of high mobility class
// ml is computed to maintain the mean. If ml is made equal to mh, both classes have same mobility
// ms is mobility spread, such that the mobility of an individual is between mi - ms and mi + ms

//static const int N = 50000;

#define N 50000
// #define num_iter 100

static const double mi = 5;  //The internal average infection rate
static const double ph = 0.8; //The fraction of high contact individuals
static const double ml = 3; //The contact rate of low risk group
static const double mh = ((mi - ml * (1 - ph)) / (ph)); //The contact rate of high risk group
static const double ms = 2; //The dispersion of contact rate
static const double T = 15000; //The total testing capacity
static const double f = 0.8; // the deflation factor for frequency of testing for high risk group
static const double F = N / T; //The average frequency of testing
static const double Fh = f * (double)N / (double)T; //The frequency of testing high risk group
static const double Fl = (F - Fh * ph) / (1-ph); //The average frequency of testing low risk group
static const double b0 = 0.025; //The base level infectivity
static const double me = 2; //External contact rate
static const double rho = 0.043; //External positivity of the environment
static const double be = rho * me * b0; //External infection rate
static double daily_tests_remaining = T; //Tests remaining for a day, will be set to T at the beginning of every cycle.
static const double eta_0 = 0.9; //efficiency of contact tracing
static const double eta_1 = 0.95; //efficiency of isolation
static const double Delay = 0; //Delay in test results
static const double sens = 0.95; //sensitivity of testing;
static const int U0 = 5; //Initial infection load
static const int Recovery_Period = 15; //Recovery Period
static const int Tmax = 120; //Maximum time period for which the program needs to run
static const bool adaptive = true;
static double tests_allocated = T;
static double last_positivity = U0/N;

// pthread_barrier_t barrier;
// pthread_mutex_t lock;
// pthread_mutex_t lock1;
//Structure definition of Agents
typedef struct agent {
	int ID;
	bool S; //Susceptible - individual status
	bool U; //Undetected - individual status
	bool P; //Positive - individual status
	bool R; //Recovered - individual status
	bool C; //In contact list
	bool T; //Tested
	int risk; //Infection - mobility risk class
	int max_contact; // Maximum number of contacts that an agent can make.
	int daily_contact; // Daily contacts - Random with a distribution
	int time_in_status; //Time spent in each status
	int time_since_test; //Time since tested
	//agent * contact;
}agent;


// Simulation structure
typedef struct Simulation {
	agent * Agents[N]; //An array of agents equal to the population.
	int Susceptible; //Total number of susceptibles
	int Undetected; //Total number of undetected infected
	int Positive; //Total number of detected infected
	int Recovered; // Total recovered
	int NewInfections; // New infections every day
}Simulation;

Simulation * Sim;

//helper function for randomizing agent indexes for contacts
int randBetween(int lower, int upper) {
	return rand() % (upper - lower + 1) + lower;
}


//Create a daily contact list as a linked list
typedef struct contact_edge{
	int id1; //Agent node 1
	int id2; //Agent node 2
	struct contact_edge * next; //Next edge - to maintain the edge path
}contact;

//contact list head and tail
contact * contact_head = NULL;
contact * contact_tail = NULL;


//Functions

//Agent functions
agent * newAgent(int id); // Constructor
void destroy_agent(agent * a); //Destructor
void increase_agent_time(int id); //Helper function for agent

//Initiatization of simulation
Simulation * initiate_simulation();

//Contact edge related finctions - Agents are the nodes
contact * create_contact_edge(int i1, int i2);
void destroy_contact_edge(contact * this);
void destroy_contact_list();
void reveal_contacts(int id);
void insert_contact(int id1, int id2);
void create_contact(int id1, int id2);
void contact_event();

//Testing functions
void testing_event(int id);
void reveal_test_results(int id);

//Initiate a day
void initialize_infections();
void initialize_day();

//Run one cycle
//1. Contact tracing testing
//2. Random testing
//3. Contact generation
void run_day();

//Recovery
void recovery(int id);
