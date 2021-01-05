#include "AgentSimulation.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

//Create a new agent with the agent ID
//Max contact is randomly drawn from a Poisson distribution with mean rate of mh for 0 and ml for 1 risk class.
//risk class is randomly assigned with a probability of ph

//Agent constructor function
agent * newAgent(int id) {
	agent * a = malloc(sizeof(agent));
	a->ID = id;
	a->S = true;
	a->U = false;
	a->P = false;
	a->R = false;
	a->C = false;
	a->T = false;
	if(((double)rand()) / ((double)RAND_MAX) > ph){
		a->risk = 0;
		a->max_contact = randBetween((int)ml - ms, (int)ml + ms);
	}else{
		a->risk = 1;
		a->max_contact = randBetween((int)mh - ms, (int)mh + ms);
	}
	a->daily_contact = 0;
	a->time_in_status = 0;
	a->time_since_test = 0;
	return a;
}

//destructor function for agents
void destroy_agent(agent * a){
	free(a);
	a = NULL;
}

//Constructor for the simulation structure
Simulation * initiate_simulation(){
	Simulation * newSimulation = malloc(sizeof(Simulation));
	newSimulation->Susceptible = N;
	newSimulation->Undetected = 0;
	newSimulation->Positive = 0;
	newSimulation->Recovered = 0;
	newSimulation->NewInfections = 0;
	time_t t;
	srand((unsigned) time(&t));
	for(int i = 0; i < N; i++){
		newSimulation->Agents[i] = newAgent(i);
	}
	return newSimulation;
}

//Destructor for the simulation structure
void delete_simulation(){
	for(int i = 0; i < N; i++){
		free(Sim->Agents[i]);
	}
	free(Sim);
}

//Create a contact edge with two agents and next contact
//Contructor for the edge of the agent graph
contact * create_contact_edge(int i1, int i2){
	contact * thisContact = malloc(sizeof(contact));
	thisContact->id1 = i1;
	thisContact->id2 = i2;
	thisContact->next = NULL;
}

//destroy a contact edge
//destructor for the edge of the agent graph
void destroy_contact_edge(contact * this){
	free(this);
	this = NULL;
}

//destroy the contact list
//Destructor for the entire graph of contacts
void destroy_contact_list(){
	contact * this = contact_head;
	while(this != NULL){
		contact * temp = this->next;
		destroy_contact_edge(this);
		this = temp;
		temp = NULL;
	}
	contact_head = NULL;
	contact_tail = NULL;
}

//Traverse the contact list and reveal contacts
//of a specific agent with ID = id
//Edge traversal for the agent graph
void reveal_contacts(int id){
	contact * this = contact_head;
	while(this != NULL){
		if(this->id1 == id){
			double r1 = ((double)rand()) / ((double)RAND_MAX);
			if(r1 < eta_0){
				Sim->Agents[this->id2]->C = true;
			}
		}
		if(this->id2 == id){
			double r1 = ((double)rand()) / ((double)RAND_MAX);
			if(r1 < eta_0){
				Sim->Agents[this->id1]->C = true;
			}
		}
		this = this->next;
	}
}

//Insert a pair of contacts in contact list
//Edge insertion
void insert_contact(int id1, int id2){
	contact * this = create_contact_edge(id1, id2);
	if(contact_head == NULL){  //contact_head == NULL
		contact_head = this;
		contact_tail = this;
	}else{
		contact_tail->next = this;
		contact_tail = this;
	}
}

//void create_contact between two individuals
//Infection transmits with probability b0
//Only if one is susceptibe and the other is undetected
//The daily contact of both individuals are updated
//The daily contact of
//Edge full constructor from agent IDs
void create_contact(int id1, int id2){
	if((Sim->Agents[id1]->daily_contact >= Sim->Agents[id1]->max_contact) ||
		(Sim->Agents[id2]->daily_contact >= Sim->Agents[id2]->max_contact) || (id1 == id2)){
			return;
	}
	if((Sim->Agents[id1]->P) || (Sim->Agents[id2]->P)){
		return;
	}
	if(((Sim->Agents[id1]->U) && ((Sim->Agents[id2]->U) || (Sim->Agents[id2]->S)))
				|| ((Sim->Agents[id2]->U) && ((Sim->Agents[id1]->U) || (Sim->Agents[id1]->S)))){
		insert_contact(id1, id2);
	}
	// if((Sim->Agents[id1]->U) && (Sim->Agents[id2]->U)){
	// 	return;
	// }
	// if((Sim->Agents[id1]->R) || (Sim->Agents[id2]->R)){
	// 	return;
	// }
	if((Sim->Agents[id1]->U) && (Sim->Agents[id2]->S)){
			if(((double)rand()) / ((double)RAND_MAX) <= b0){
				Sim->Agents[id2]->S = false;
				Sim->Agents[id2]->U = true;
				Sim->Agents[id2]->P = false;
				Sim->Agents[id2]->R = false;
				Sim->Agents[id2]->time_in_status = 0;
				(Sim->Susceptible)--;
				(Sim->Undetected)++;
				(Sim->NewInfections)++;
			}
	}
	if((Sim->Agents[id2]->U) && (Sim->Agents[id1]->S)){
			if(((double)rand()) / ((double)RAND_MAX) < b0){
				Sim->Agents[id1]->S = false;
				Sim->Agents[id1]->U = true;
				Sim->Agents[id1]->P = false;
				Sim->Agents[id1]->R = false;
				Sim->Agents[id1]->time_in_status = 0;
				(Sim->Susceptible)--;
				(Sim->Undetected)++;
				(Sim->NewInfections)++;
			}
	}
	(Sim->Agents[id1]->daily_contact)++;
	(Sim->Agents[id2]->daily_contact)++;
}

//Helper function for the agents
void increase_agent_time(int id){
	(Sim->Agents[id]->time_in_status)++;
	(Sim->Agents[id]->time_since_test)++;
}



//Draw two random ids of agents
//Generate contact of the two agents
//Created full contact event
void contact_event(){
	int id1 = randBetween(0, N - 1);
	if(id1 < 0){
		id1 = 0;
	}
	if(id1 >= N){
		id1 = N - 1;
	}
	int id2 = randBetween(0, N - 1);
	if(id2 < 0){
		id2 = 0;
	}
	if(id2 >= N){
		id2 = N - 1;
	}
	create_contact(id1, id2);
}

//Testing event
//We create a test for an individual who is not tested within the past
//max frequency number of days
//reduce daily test capacity by 1 if tested
//Testing for all but positive cases with known frequency
void testing_event(int id){
	if(Sim->Agents[id]->P){
		return;
	}
	if(daily_tests_remaining <= 0){
		return;
	}
	if(!(Sim->Agents[id]->T) && (Sim->Agents[id]->time_since_test >= Fh)
							&& (Sim->Agents[id]->risk == 1)){
		Sim->Agents[id]->T = true;
		Sim->Agents[id]->time_since_test = 0;
		daily_tests_remaining--;
		return;
	}
	if(!(Sim->Agents[id]->T) && (Sim->Agents[id]->time_since_test >= Fl)
							&& (Sim->Agents[id]->risk == 0)){
		Sim->Agents[id]->T = true;
		Sim->Agents[id]->time_since_test = 0;
		daily_tests_remaining--;
		return;
	}
}


//Reveal test results if the testing is true and
//the number of days since testing is >= delay in results
void reveal_test_results(int id){
	if(!(Sim->Agents[id]->U)){
		Sim->Agents[id]->T = false;
		return;
	}
	if(!(Sim->Agents[id]->T)){
		return;
	}
	if((Sim->Agents[id]->time_since_test < Delay)){
		return;
	}
	double r1 = ((double)rand()) / ((double)RAND_MAX);
	double r2 = ((double)rand()) / ((double)RAND_MAX);
	Sim->Agents[id]->T = false;
	if((r1 <= sens) && (r2 <= eta_1)){
		Sim->Agents[id]->U = false;
		Sim->Agents[id]->P = true;
		Sim->Agents[id]->time_in_status = 0;
		(Sim->Undetected)--;
		(Sim->Positive)++;
		reveal_contacts(id);
	}
}

//Destroy contact list
//Initialize testing capacity
//Resets the agent contact list
void initialize_day(){
	if(adaptive){
		double temp = ((double)Sim->NewInfections) / tests_allocated;
		tests_allocated = (1.3 * temp / (last_positivity + 0.001)) * tests_allocated;
		last_positivity = temp;
		// fprintf(stderr, "%f : %f\n", tests_allocated, last_positivity);
		if(tests_allocated > T){
			tests_allocated = T;
		}
		if(tests_allocated < 500){
			tests_allocated = 500;
		}
		daily_tests_remaining = tests_allocated;
	}else{
		daily_tests_remaining = T;
	}
	Sim->NewInfections = 0;
	destroy_contact_list();
	for(int i = 0; i < N; i++){
		Sim->Agents[i]->daily_contact = 0;
	}
}


//Initialize the number of infections
//Randomly initialize the first U0 number of infections
void initialize_infections(){
	for(int i = 0; i < U0; i++){
		int id = randBetween(0, N-1);
		if(id < 0){
			id = 0;
		}
		if(id >= N){
			id = N - 1;
		}
		Sim->Agents[id]->S = false;
		Sim->Agents[id]->U = true;
		Sim->Agents[id]->P = false;
		Sim->Agents[id]->R = false;
		(Sim->Susceptible) -= 1;
		(Sim->Undetected) += 1;
	}
}

//Recover an infected individual after Recovery Period
//Transition function from infected to recovery
void recovery(int id){
	double thisRecovery = 1/((double)randBetween(Recovery_Period - 4, Recovery_Period + 4));
	double r = ((double)rand()) / ((double)RAND_MAX);
	if(Sim->Agents[id]->U){
		if(thisRecovery >= r){
			Sim->Agents[id]->S = false;
			Sim->Agents[id]->U = false;
			Sim->Agents[id]->P = false;
			Sim->Agents[id]->R = true;
			(Sim->Recovered)++;
			(Sim->Undetected)--;
		}
	}
	if(Sim->Agents[id]->P){
		if(thisRecovery >= r){
			Sim->Agents[id]->S = false;
			Sim->Agents[id]->U = false;
			Sim->Agents[id]->P = false;
			Sim->Agents[id]->R = true;
			(Sim->Recovered)++;
			(Sim->Positive)--;
		}
	}
	// if(Sim->Agents[id]->U){
	// 	if(Sim->Agents[id]->time_in_status >= Recovery_Period){
	// 		Sim->Agents[id]->U = false;
	// 		Sim->Agents[id]->R = true;
	// 		(Sim->Recovered)++;
	// 		(Sim->Undetected)--;
	// 	}
	// }
	// if(Sim->Agents[id]->P){
	// 	if(Sim->Agents[id]->time_in_status >= Recovery_Period){
	// 		Sim->Agents[id]->P = false;
	// 		Sim->Agents[id]->R = true;
	// 		(Sim->Recovered)++;
	// 		(Sim->Positive)--;
	// 	}
	// }
}

//External infection
//Based on the external positivity
void external_infection(int id){
	if(Sim->Agents[id]->S){
		double r1 = ((double)rand()) / ((double)RAND_MAX);
		if(r1 < be){
			Sim->Agents[id]->S = false;
			Sim->Agents[id]->U = true;
			(Sim->Undetected)++;
			(Sim->Susceptible)--;
		}
	}
}

//Run a day of simulation
//Initialize day
//Perform contact testing
//Perform testing
//Perform random contacts
//Reveal contacts
void run_day(){
	initialize_day();
	//Contact Testing
	for(int i = 0; i < N; i++){
		if(Sim->Agents[i]->C){
			// fprintf(stderr, "Contact Traced ID %d\n", Sim->Agents[i]->ID);
			testing_event(i);
			Sim->Agents[i]->C = false;
		}
	}
	for(int i = 0; i < N; i++){
		// int id = randBetween(0, N-1);
		// if(id < 0){
		// 	id = 0;
		// }
		// if(id >= N){
		// 	id = N - 1;
		// }
		testing_event(i);
		if(daily_tests_remaining <= 0){
			break;
		}
	}
	// fprintf(stderr, "I AM HERE\n");
	// for(int i = 0; i < N; i++){
	// 	fprintf(stderr, "STARTING AGENT %d\n", i);
	// 	while((Sim->Agents[i]->daily_contact) < (Sim->Agents[i]->max_contact)){
	// 		int id = randBetween(0, N-1);
	// 		fprintf(stderr, "CONTACT ID - %d \nDAILY CONTACTS FOR AGENT %d AND MAX CONTACT IS %d\n", id, Sim->Agents[i]->daily_contact, Sim->Agents[i]->max_contact);
	// 		create_contact(i, id);
	// 	}
	// }
	for(int i = 0; i < mi * N; i++){
		contact_event();
	}
	for(int i = 0; i < N; i++){
		external_infection(i);
	}
	for(int i = 0; i < N; i++){
		reveal_test_results(i);
		increase_agent_time(i);
		recovery(i);
	}
}


// void * run_all_days(void * fname){
// 	FILE * file = NULL;
// 	// int num = 100;
//  // pthread_barrier_wait(&barrier);
// 	// pthread_mutex_lock(&lock);
// 	file = fopen((char *)fname, "w");
// 	// pthread_mutex_unlock(&lock);
// 	fprintf(file, "Susceptibles,Undetected,Positive,Recovered\n");
// 	Sim = initiate_simulation();
// 	initialize_infections();
// 	for(int i = 0; i < Tmax; i++){
// 		run_day();
// 		if(file){
// 			fprintf(file, "%d,%d,%d,%d\n",
// 				Sim->Susceptible, Sim->Undetected, Sim->Positive, Sim->Recovered);
// 		}else{
// 			fprintf(stdout, "Susceptibles: %d; Undetected: %d; Positive: %d; Recovered: %d\n",
// 				Sim->Susceptible, Sim->Undetected, Sim->Positive, Sim->Recovered);
// 		}
// 	}
// 	delete_simulation();
// 	if(file){
// 		fclose(file);
// 	}
// 	// pthread_barrier_wait(&barrier);
// 	return NULL;
// }


int main(int argc, char ** argv) {
	if(argc <= 1){
		return 0;
	}
	FILE * file = NULL;
	int num = 100;
	// if(argc > 1){
	// 	fprintf(stderr, "%d : %s : %s\n", argc, argv[0], argv[1]);
	// 	file = fopen(argv[1], "w");
	// 	if(file == NULL){
	// 		fprintf(stderr, "ERROR OPENING FILE\n");
	// 	}
	// 	// if(argc == 3){
	// 	// 	num = atoi(argv[2]);
	// 	// }
	// }
	// clock_t begin = clock();

	char * fnames;
	// pthread_t threads[num_iter];

	// pthread_barrier_init(&barrier, NULL, num_iter+1);

	while(num > 0){
		fprintf(stderr, "STARTING NEXT ROUND %d\n", num);
		// pthread_mutex_lock(&lock);
		if(argc > 1){
			fnames = calloc(1000, 1);
			sprintf(fnames,"./SimulationOutput/%d_HEADER_%s",num,argv[1]);
			puts(fnames);
			file = fopen(fnames, "w");
			memset(fnames,'\0',1000);

		if(file){
			fprintf(file, "\n SIMULATION PARAMETERS \n");
			fprintf(file, "Total Population : %d\n", N);
			fprintf(file, "Total Tests : %f \n", T);
			fprintf(file, "Average Internal Mobility : %f \n", mi);
			fprintf(file, "Average External Mobility : %f \n", me);
			fprintf(file, "Base Infectivity : %f \n", b0);
			fprintf(file, "External Positivity : %f \n", rho);
			fprintf(file, "Isolation Efficiency : %f \n", eta_1);
			fprintf(file, "Contact Tracing Efficiency : %f \n", eta_0);
			fprintf(file, "Testing Delay : %f \n", Delay);
			fprintf(file, "Test Sensitivity : %f \n", sens);
		}else{
			fprintf(stderr, "\n SIMULATION PARAMETERS \n");
			fprintf(stderr, "Total Population : %d\n", N);
			fprintf(stderr, "Total Tests : %f \n", T);
			fprintf(stderr, "Average Internal Mobility : %f \n", mi);
			fprintf(stderr, "Average External Mobility : %f \n", me);
			fprintf(stderr, "Base Infectivity : %f \n", b0);
			fprintf(stderr, "External Positivity : %f \n", rho);
			fprintf(stderr, "Isolation Efficiency : %f \n", eta_1);
			fprintf(stderr, "Contact Tracing Efficiency : %f \n", eta_0);
			fprintf(stderr, "Testing Delay : %f \n", Delay);
			fprintf(stderr, "Test Sensitivity : %f \n", sens);
		}
		fclose(file);
		sprintf(fnames,"./SimulationOutput/%d_OUTPUT_%s",num,argv[1]);
		puts(fnames);
		file = fopen(fnames, "w");
		free(fnames);
		fprintf(file, "Susceptibles,Undetected,Positive,Recovered\n");
	}
		Sim = initiate_simulation();
		initialize_infections();
		for(int i = 0; i < Tmax; i++){
			run_day();
			if(file){
				fprintf(file, "%d,%d,%d,%d, %f\n",
					Sim->Susceptible, Sim->Undetected, Sim->Positive, Sim->Recovered, tests_allocated);
			}else{
				fprintf(stdout, "Susceptibles: %d; Undetected: %d; Positive: %d; Recovered: %d\n",
					Sim->Susceptible, Sim->Undetected, Sim->Positive, Sim->Recovered);
			}
		}
		delete_simulation();
		// pthread_create(&threads[num - 1], NULL, run_all_days, (void *)fnames[num - 1]);
		// pthread_mutex_unlock(&lock);
		// run_all_days((void *)fnames[num - 1]);

		num--;
		if(file){
			fclose(file);
		}
	}

	// pthread_barrier_wait(&barrier);
	// pthread_barrier_wait(&barrier);


	// for(int i = 0; i < num_iter; i++){
	// 	// pthread_join(threads[i], NULL);
	// 	free(fnames[i]);
	// 	fnames[i] = NULL;
	// }
	// clock_t end = clock();
	// double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	// fprintf(stderr, "%f\n", time_spent);
	// int Sus = 0;
	// int Infec = 0;
	// int Pos = 0;
	// int rec = 0;
	// for(int i = 0; i < N; i++){
	// 	if(Sim->Agents[i]->S)Sus++;
	// 	if(Sim->Agents[i]->U)Infec++;
	// 	if(Sim->Agents[i]->P)Pos++;
	// 	if(Sim->Agents[i]->R)rec++;
	// }
	// fprintf(stderr, "Susceptible : %d | Undetected : %d | Positive : %d | Recovered : %d\n", Sus, Infec, Pos, rec);
	// pthread_barrier_destroy(&barrier);
	return 0;
}
