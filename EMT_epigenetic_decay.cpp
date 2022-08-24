////////////////////////////////////////////////////////////////////////////////////////////////
////                          	  EMT Epigenetic decay                                     /////
////                        Jain et al., Date: 23th Aug 2022                              /////
//////////////////////////////////////////////////////////////////////////////////////////////// 
#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <string>
#include <fstream>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/numeric/odeint.hpp>

#define NUMNODES 6

typedef boost::array<double, NUMNODES> cell_state; // a cells' variable vector for ODE simulation function {miR-200, mRNA ZEB, ZEB, SNAIL, Z0u200, S0}
typedef boost::unordered_map<int, boost::array<double, NUMNODES>> population; // stores the values of 6 variables for every cell in the population


std::mt19937 generator; // for random number generation

double temp_var[10]; // temp_var hold the value for initialization parameters provided at start of simulation

// defining global variables - beta: time constant for z0u200 threshold changes; gamma_SNAIL: SNAIL's return rate; alpha: epigenetic regulation strength
double beta, alpha, gamma_SNAIL;


//////////////////////// k combination out of n //////////////////////////////
double nchoosek(int n, int k)
{
	double result = 1.0;
	double result0 = 1.0;

	for(int i = 0; i < k; i++)
	{
		result *= (n - i);
	}
	for(int i = 1; i <= k; i++)
	{
		result0 *= i;
	}

	return(result / result0);
}
///////////////////////////////////////////////////////////////////////////////////

// function to break a string and separate numeral character (initializers) based on the delimiter
void tokenize(std::string s, std::string del = " ")
{
	int start = 0;
	int end = s.find(del);
	int i = 0;
	while (end != -1) {
		temp_var[i] = std::stod(s.substr(start, end - start));
		start = end + del.size();
		end = s.find(del, start);
		i++;
	}
	temp_var[i] = std::stod(s.substr(start, end - start));
}
///////////////////////////////////////////////////////////////////////////////////

// Function to integrate EMT network differential equation
// takes input as reference of cell population vector array - P, standard deviation of noise in SNAIL's SDE - sigma, and euler integration step size (in units of hours) : euler_step
// output is updated population vector array - P - at the next time step (t + euler_step)
void euler_integration_step(population & P, double sigma, double euler_step) 

{

	std::normal_distribution <> dist{ 0.0, 1.0 }; // initialization normal random distribution with mean 0 and standard deviation 1 

	double temp_variable[6]; // temporary variable to hold cell state variables in the sequence {miR-200, mRNA ZEB, ZEB, SNAIL, Z0u200, S0}

	// initializing system's ODE parameters

	// Degradation rates
	double ku200 = 0.05, kmz = 0.5, kz = 0.1;
	 
	// Transcription rate:
	double gu200 = 2100, gmz = 11, gz = 100;

	// Hills function threshold :
	double z0u200_0 = 220000, z0mz = 25000, s0u200 = 180000, s0mz = 180000, u2000 = 10000;

	// Cooperativity:
	double nzu200 = 3, nsu200 = 2, nzmz = 2, nsmz = 2, nu200 = 6;

	// fold change
	double lamdazu200 = 0.1, lamdasu200 = 0.1, lamdazmz = 7.5, lamdasmz = 10;

	// Euler intregration of ODE
	for(int i = 0; i < P.size(); i++){

		double Mu0=1/std::pow((1+P[i][0]/u2000),nu200);
		double Mu1=std::pow((P[i][0]/u2000),1)/std::pow((1+P[i][0]/u2000),nu200);
		double Mu2=std::pow((P[i][0]/u2000),2)/std::pow((1+P[i][0]/u2000),nu200);
		double Mu3=std::pow((P[i][0]/u2000),3)/std::pow((1+P[i][0]/u2000),nu200);
		double Mu4=std::pow((P[i][0]/u2000),4)/std::pow((1+P[i][0]/u2000),nu200);
		double Mu5=std::pow((P[i][0]/u2000),5)/std::pow((1+P[i][0]/u2000),nu200);
		double Mu6=std::pow((P[i][0]/u2000),6)/std::pow((1+P[i][0]/u2000),nu200);


			
		double Hillszu200=(1+lamdazu200*std::pow((P[i][2]/P[i][4]),nzu200))/(1+std::pow((P[i][2]/P[i][4]),nzu200));
		double Hillssu200=(1+lamdasu200*std::pow((P[i][3]/s0u200),nsu200))/(1+std::pow((P[i][3]/s0u200),nsu200));
		double Hillszmz=(1+lamdazmz*std::pow((P[i][2]/z0mz),nzmz))/(1+std::pow((P[i][2]/z0mz),nzmz));
		double Hillssmz=(1+lamdasmz*std::pow((P[i][3]/s0mz),nsmz))/(1+std::pow((P[i][3]/s0mz),nsmz));

		  

		temp_variable[0] = P[i][0] + (gu200*Hillszu200*Hillssu200-P[i][1]*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-ku200*P[i][0])* euler_step;
		temp_variable[1] = P[i][1] + (gmz*Hillszmz*Hillssmz-P[i][1]*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmz*P[i][1]) * euler_step;
		temp_variable[2] = P[i][2] + (gz*P[i][1]*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-kz*P[i][2]) * euler_step;
		temp_variable[4] = P[i][4] + ((1/beta) * (z0u200_0 - P[i][4] - alpha*P[i][2])) * euler_step ; 
		temp_variable[5] = P[i][5];



		while(true){

			// Euler integration of Stochastic Differential Equation
			temp_variable[3] = P[i][3] + (gamma_SNAIL*(1 - P[i][3]/P[i][5]))* euler_step + dist(generator) * sigma * std::pow(euler_step,0.5); 

			if(temp_variable[3] > 0)
				break;
		}

		P[i][0] = temp_variable[0]; 
		P[i][1] = temp_variable[1];
		P[i][2] = temp_variable[2];
		P[i][3] = temp_variable[3];
		P[i][4] = temp_variable[4];
		P[i][5] = temp_variable[5];

		
	}
}


/////////////////////////// to assign phenotypes to cells /////////////////////////
/////////////////////////// based upon levels of mZEB levels //////////////////////
// takes input as reference of population vector array (P) and its corresponding cell phenotype array mapping
// output is updated phenotype of the population
void get_phenotypes(population &P, boost::unordered_map <int, int> &phenotype)
{
	int state = -1;

	for(int i =0; i < P.size(); i++){
		
		if(P[i][1] < 160){
			
			state = 0;
		}
		else if (P[i][1] >= 160 && P[i][1] <= 568)
		{
			state = 1;
		} 
		else{
			state = 2;
		}


		if(state == -1)
		{
			std::cout << "Error in phenotype assignment." << "\n";
		}
		else
		{
			phenotype[i] = state;
		}
	}

}

// function either initialize cell variables at the beginning of simulation, sets S0 levels to S02 at the start of EMT induction or sets S0 to S01 at the end of EMT induction period
// inputs are reference to the population vector array (P), total cell population size (pop_size), pre-induction SNAIL level (S01), pre-induction SNAIL level (S02),...
// time constant for accumulation of epigenetic changes during EMT (beta_for), time constant for decay of epigenetic changes during withdrawal (beta_rev), SNAIL's return rate during EMT induction (gamma_SNAIL_for),...
// SNAIL's return rate during withdrawal (gamma_SNAIL_rev), euler integration step size (euler_step), standard deviation of noise in SNAIL's SDE (sigma)

// sampling_number = 1 denotes that simulation is started and cell variables has to be initialized
// sampling_number = 2 denotes that EMT induction is starting and therefore, the function sets et S0 levels to S02
// sampling_number = 3 denotes that period of EMT induction has ended and withdrawal is started. Therefore, set S0 levels to S01

void SNAIL_normal(population& P, const int pop_size, int sampling_number, int S01,  int S02, double beta_for, double beta_rev, double gamma_SNAIL_for, double gamma_SNAIL_rev, double euler_step, double sigma){

	std::normal_distribution <> dist{ 0.0, 1.0 };

	int count = 0;
	double temp_SNAIL;

	while(count < pop_size){
		
		if(sampling_number == 1){
			
			
			while(true){

				temp_SNAIL = S01 + sigma*std::pow((S01/(2*gamma_SNAIL_rev)),0.5) * dist(generator);

				if(temp_SNAIL > 0)
					break;
			}

			P[count][3] = temp_SNAIL;

			P[count][5] = temp_SNAIL; // this will maintain the SNAIL distribution characteristic sampled above during cells initialization
			beta = beta_rev;
			gamma_SNAIL = gamma_SNAIL_rev;
		}
		
		else if(sampling_number == 2)
		{
			P[count][5] = S02;	
			beta = beta_for;
			gamma_SNAIL = gamma_SNAIL_for;
		}
		else
		{
			P[count][5] = S01;
			beta = beta_rev;
			gamma_SNAIL = gamma_SNAIL_rev;
		}

		count++;
	}

	count = 0;
	if(sampling_number == 1){
		
		while(count < pop_size)

		{
		

			P[count][0] = 0.0;
			P[count][1] = 0.0;
			P[count][2] = 0.0;
			P[count][4] = 220000;

			count++;

		}

		double t = 0.0;
		
		// to stabilize all the variable value according to the sampled SNAIL level
		while(t < 1000.0){
			euler_integration_step(P,0.0,euler_step); // here, the noise level are zero
			t = t + euler_step; 
		}
		///////////////////////////////////////////////////////////////////////////

		// initial burn in time with fluctuation of cellular SNAIL level around the population mean S0 
		
		for(int i = 0; i < P.size(); i++){
			P[i][5] = S01;
		}

		t = 0.0;
		
		while(t < 1000.0){
			euler_integration_step(P,sigma,euler_step);
			t = t + euler_step; 
		}
		///////////////////////////////////////////////////////////////////////////////////////////////
	}

}

////////////////////////// Simulating EMT induction and withdrawal ////////////////////////////////////////////////////////////
// inputs are reference to the population vector array (P), total cell population size (pop_size), pre-induction (post-withdrawal) population mean SNAIL level (S01), pre-induction population mean SNAIL level (S02),...
// time constant for accumulation of epigenetic changes during EMT (beta_for), time constant for decay of epigenetic changes during withdrawal (beta_rev), SNAIL's return rate during EMT induction (gamma_SNAIL_for),...
// SNAIL's return rate during withdrawal (gamma_SNAIL_rev), euler integration step size (euler_step), standard deviation of noise in SNAIL's SDE (sigma)...
// total simulation time in hours (end_time), iteration number of the simulation (sim_num), pointer to the variable treat_period that contains EMT induction and withdrawal duration...
// CV of SNAIL distribution (CV), S0 values for which the above CV is defined (S0_CV), 50% decorrelation time of stochastic SNAIL trajectory (tau)

// function stores the phenotype of the cells in the population, and their corresponding cell state variables at the end of days. All data are stored in .csv file format.

int simulate_normal(population& P, double end_time, int pop_size, int sim_num, int * treat_period, double beta_for,double beta_rev , double gamma_SNAIL_for, double gamma_SNAIL_rev, int S01, int S02, double sigma, double euler_step, double S0_CV, double CV, double tau)
{

	boost::unordered_map <int, int> phenotype;
	boost::array<double, 3> count;

	double day_indx = 1;

	double t = 0.0; // time is in hours
	int flag = 0; // flag variable is used to switch S0 to S01 at the end of EMT induction

	std::ofstream end_of_day, SNAIL_distribution, mZEB_distribution, Z0u200_distribution , S0_distribution, miR200_distribution, ZEB_distribution;

	// population phenotype distribution data of multiple simulation are appended in the same file. 
	//However, each cell state variable data is saved in distinct file for every new simulation run

	end_of_day.open("Data/end_of_day_counts_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev)) + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau))  + ".csv", std::ios::app);
	
	miR200_distribution.open("Data/miR200_dist_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + "_" + std::to_string(sim_num) + ".csv");

	mZEB_distribution.open("Data/mZEB_dist_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + "_" + std::to_string(sim_num) + ".csv");
	
	ZEB_distribution.open("Data/ZEB_dist_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + "_" + std::to_string(sim_num) + ".csv");

	SNAIL_distribution.open("Data/snail_dist_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + "_" + std::to_string(sim_num) + ".csv");

	Z0u200_distribution .open("Data/threshold_dist_"  + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + "_" + std::to_string(sim_num) + ".csv");

	S0_distribution.open("Data/final_snail_dist_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + "_" + std::to_string(sim_num) + ".csv");
	
	get_phenotypes(P, phenotype);

	count = {0.0, 0.0, 0.0};
	
	for(int i = 0; i < P.size(); i++)
	{
		count[phenotype[i]] += 1.0;
	}
	
	//// the phenotypic distribution is stored each day in absolute cell counts ////
	end_of_day << "0" << "," << count[0] << "," << count[1] << "," << count[2] << "," << P.size() << "\n";

	for(int i = 0; i < P.size(); i++){
		
		if(i < P.size() -1){

			miR200_distribution << P[i][0] << ",";
			mZEB_distribution << P[i][1] << ",";
			ZEB_distribution << P[i][2] << ",";
			SNAIL_distribution << P[i][3] << ",";
			Z0u200_distribution  << P[i][4] << ",";
			S0_distribution << P[i][5] << ",";			
		}
		else{
			miR200_distribution << P[i][0] << "\n";
			mZEB_distribution << P[i][1] << "\n";
			ZEB_distribution << P[i][2] << "\n";
			SNAIL_distribution << P[i][3] << "\n";
			Z0u200_distribution  << P[i][4] << "\n";
			S0_distribution << P[i][5] << "\n";	
		}	

	}

	
	gamma_SNAIL = gamma_SNAIL_for;
	beta = beta_for;


	// set S0 to S02 to induce EMT in cells
	SNAIL_normal(P,pop_size,2, S01, S02,beta_for, beta_rev, gamma_SNAIL_for, gamma_SNAIL_rev,  euler_step, sigma);


	while(t <= end_time)
	{

		if(t > (24 * treat_period[0]) && flag == 0){
			// EMT induction period is over and withdrawal is started
			gamma_SNAIL = gamma_SNAIL_rev;
			beta = beta_rev;

			// set S0 to S02 to induce EMT in cells
			SNAIL_normal(P,pop_size,3, S01, S02,beta_for,beta_rev, gamma_SNAIL_for, gamma_SNAIL_rev, euler_step, sigma);
			flag = 1;
		}	
		
		euler_integration_step(P,sigma, euler_step); // increment entire population dynamics by one euler step (defined by euler_step), sigma is the noise level in SNAIL dynamics 
		
		t = t + euler_step; // here, time is in hours


		if(true)
		{

			if(t >= 24* day_indx){

				get_phenotypes(P, phenotype);

				count = {0.0, 0.0, 0.0};
				
				for(int i = 0; i < P.size(); i++)
				{
					count[phenotype[i]] += 1.0;
				}
				
				//// the phenotypic distribution is stored each day in absolute cell counts ////
				end_of_day << day_indx << "," << count[0] << "," << count[1] << "," << count[2] << "," << P.size() << "\n";
				

				for(int i = 0; i < P.size(); i++){
					
					if(i < P.size() -1){
						miR200_distribution << P[i][0] << ",";
						mZEB_distribution << P[i][1] << ",";
						ZEB_distribution << P[i][2] << ",";
						SNAIL_distribution << P[i][3] << ",";
						Z0u200_distribution  << P[i][4] << ",";				
						S0_distribution << P[i][5] << ",";
					}
					else{
						miR200_distribution << P[i][0] << "\n";
						mZEB_distribution << P[i][1] << "\n";
						ZEB_distribution << P[i][2] << "\n";
						SNAIL_distribution << P[i][3] << "\n";
						Z0u200_distribution  << P[i][4] << "\n";				
						S0_distribution << P[i][5] << "\n";
					}	

				}

				day_indx++;

			}
		}
	}

	if(true)
	{
		end_of_day.close();
		miR200_distribution.close();
		mZEB_distribution.close();
		ZEB_distribution.close();
		SNAIL_distribution.close();	
		Z0u200_distribution .close();
		S0_distribution.close();
	
	}

	return(0);
}



int main(int argc, char* argv[])
{
	 
	tokenize(argv[1], "_"); /// to separate parameters from input string 

	boost::unordered_map <int, int> phenotype; // cell population' phenotype vector array
	double euler_step = 0.1; // euler integration step size

	// take input CV, and at what S0 it corresponds to; tau(50 percent decoorelation time)
	// input to be given while running the code is in the following sequence

	// induction_days_withdrawal_days_beta_for_beta_rev_alpha_S0_CV_CV_S01_S02_tau
	// e.g. 39_100_240_720_0.15_150000_0.3_100000_350000_120

	int treat_period[2] = {temp_var[0], temp_var[1]}; // induction period, withdrawal period
	double beta_for = temp_var[2]; // time constant for accumulation of epigenetic changes during EMT (beta_for). Initialize beta_for to 1 to simulate no epigenetic regulation scenario 
	double beta_rev = temp_var[3]; // time constant for decay of epigenetic changes during withdrawal (beta_rev). Initialize beta_rev to 1 to simulate no epigenetic regulation scenario
	alpha = temp_var[4]; // epigenetic regulation strength. Initialize alpha to 0 to simulate no epigenetic regulation scenario 
	double S0_CV = temp_var[5]; // mean of the distribution having below specified CV 
	double CV = temp_var[6]; // CV of dsitribution whose mean is S0_CV
	int S01 = temp_var[7]; // pre-induction (post-withdrawal) population mean SNAIL levels 
	int S02 = temp_var[8]; // post-induction population mean SNAIL levels
	double tau = temp_var[9]; // 50 percent decorrelation time in hrs

	double gamma_SNAIL_for = std::log(2)*S02/tau; // SNAIL's return rate during EMT induction (gamma_SNAIL_for)
	double gamma_SNAIL_rev = std::log(2)*S01/tau; // SNAIL's return rate during withdrawal (gamma_SNAIL_rev)

	double sigma = CV*S0_CV*std::pow((2*std::log(2)/tau),0.5); // standard deviation of noise in SNAIL's SDE (sigma). Put sigma = 0.0 to produce determinitic results 

	int sim_num; // iteration number of the simulation
	int pop_size = 10000; // cell population size in a simulation run. Put pop_size = 1 to produce deterministic results
	int total_inde_runs = 3; // total number of simulation runs. Put total_inde_runs = 1 to produce deterministic results
	int time_in_days = treat_period[0] + treat_period[1] ; // total simulation time
	int file_rows = 0;

	std::ifstream csv_file;

	std::string line;

	while(true){
		
		// create 'Data' folder in the dicrectory where this code file is saved, if you are renning the code for the first time 
		//// check how many runs of simulation data is stored in phenotype distribution data file in the 'Data' folder 
		csv_file.open("Data/end_of_day_counts_" + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(int(beta_for)) + "_" + std::to_string(int(beta_rev))  + "_" + std::to_string(alpha) + "_" + std::to_string(int(S0_CV)) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(int(tau)) + ".csv");

		
		file_rows = 0;		
		while(getline(csv_file,line)){
			file_rows++;
		}

		csv_file.close();
		
		/////////// don't proceed if file has data of required number of runs or incomplete previous simulation data  
		if(file_rows >= (total_inde_runs * (treat_period[0]+treat_period[1]+1)) || (file_rows%(treat_period[0]+treat_period[1]+1)) != 0){
			
			if(file_rows%(treat_period[0]+treat_period[1]+1) != 0)
				std::cout << "stopped bacause of already missing data \n";
			else if(file_rows > (total_inde_runs * (treat_period[0]+treat_period[1]+1)))
				std::cout << total_inde_runs * (treat_period[0]+treat_period[1]+1) << " enough data already present\n";
			break;
		} 
		//////////////////////////////////////////////////////////////////////////////////////
		
		sim_num = file_rows/(treat_period[0]+treat_period[1]+1);

		generator = std::mt19937(std::time(NULL) + file_rows); // std::time(NULL) + file_rows is the initial seed to generate random number

		population P; // define population vector array P

		// initialize cell population with SNAIL and other cell state variables values
		SNAIL_normal(P,pop_size,1, S01, S02, beta_for, beta_rev, gamma_SNAIL_for, gamma_SNAIL_rev, euler_step, sigma); // change initial pop size (200 now) as 1 for single cell simualtion
	
		std::cout << "initialization complete for " + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(beta_for) + "_" + std::to_string(beta_rev) + "_" + std::to_string(alpha) + "_" + std::to_string(S0_CV) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) + "_" + std::to_string(tau) << std::endl;
		
		// simulaate EMT induction and inducer withdrawal
		simulate_normal(P, 24* (treat_period[0]+treat_period[1]), pop_size, sim_num, treat_period, beta_for, beta_rev , gamma_SNAIL_for, gamma_SNAIL_rev, S01, S02, sigma, euler_step, S0_CV, CV, tau);

		std::cout << "simulation complete for " + std::to_string(treat_period[0]) + "_" + std::to_string(treat_period[1]) + "_" + std::to_string(beta_for) + "_" + std::to_string(beta_rev)  
		+ "_" + std::to_string(alpha) + "_" + std::to_string(S0_CV) + "_" + std::to_string(CV) + "_" + std::to_string(S01) + "_" + std::to_string(S02) 
		<< "_" <<std::to_string(int(tau)) << "_" << std::to_string(sim_num) << std::endl;

		P.clear(); // clear the cell population
	}
	
	return(0);
	
}