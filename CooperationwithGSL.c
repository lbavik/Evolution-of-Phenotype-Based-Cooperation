/* This code was written by Linnea Bavik during her Spring 2019 research rotation at Emory University in the Weissman lab.
It describes a model for the evolution of phenotype-based cooperation strategies in a population.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>

/*These set parameters for the simulation, note not all of these are necessarily in use, for instance, PEOPLE_WATCHED refers
to the number of members of the population individuals sample to make decisions about their cooperation. Currently they simply sample
the entire population so this is an unused parameter. Similarly, PARTNERS is unused, as every individual plays with every other individual. 

The April 4 2019 version of this code refers to population evolution in an invasion scenario when individuals use either a strict or
loose hard threshold of phenotype distance to decide cooperation or defection. The game matrix gives benefit B for mutual cooperation,
benefit C for mutual defection, B+C for defecting upon one's partner and 0 for being defected on.

 --------->  Note: HARDCODEDTHRESHOLD2 refers to the larger (the generous/ strategy type +1/ cooperator's) threshold <-----
*/


int POPULATION = 1000;
int GENERATIONS = 5000;
int PARTNERS = 5;
int PHENOTYPES = 1;
double SELECTIONSTRENGTH = 1;
double PHENOTYPE_MUTATION = .5;
double STRATEGY_MUTATION = .02;

int INITIAL_POPULATION = -1;
int BURN_IN = 2000;

double COST = 1;
double BENEFIT = 1.1;

int PEOPLE_WATCHED = 10;
double ACCEPTANCE_RANGE = .8;
double HARDCODEDTHRESHOLD1 = 5;
double HARDCODEDTHRESHOLD2 = 10;

int DISCRETE_DIST_PRECISION = 4;


void interaction(double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, int *mutual_coop, int *mutual_defec, int *duped);

void newgeneration(const gsl_rng * r, double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time);


#define num_phenos 1
#define num_strats 1
#define num_benefits 1
#define num_thresholds1 1
#define num_thresholds2 1

/* These are arrays of values one wishes to test. This is very easily parallelizable, but hasn't been done yet. */


double pheno_mutations[num_phenos] = {.5};
double strat_mutations[num_strats] = {.01};
double benefits[num_benefits] = {1.5};
double smaller_thresholds[num_thresholds1] = {.1};
double larger_thresholds[num_thresholds2] = {1};


int main()
{

/*For GSL Library multinomial draw */

const gsl_rng_type *T; 
/*generator type */

gsl_rng * r; 
/*creates a random number generator */

gsl_rng_env_setup(); 
/*read from environment variable */

//gsl_rng_default = gsl_rng_mt19937;
gsl_rng_default_seed = time(NULL);

T = gsl_rng_default; 
/* choose defult generator type which I set to the mersenne twister. */

gsl_rng_default = gsl_rng_mt19937;

r = gsl_rng_alloc (T); 
/*creates an instance */

        srand(time(NULL));


        double members[POPULATION][PHENOTYPES+1];
        double results[POPULATION];



//printf("STR_MUT, PHE_MUT, B, THRES1, THRES2, time, COOPs, DEFs, COOP_VAR, DEF_VAR, TOTAL_VAR, COOP_MEAN, DEF_MEAN, TOTAL_MEAN, CC_GAMES, DD_GAMES, CD_GAMES\n");



/*HERE's THE BIG OL LOOP */
/* (Please forgive the atrociously-named for loop variables) */

for (int xx = 0; xx < num_phenos; xx++)
{

        PHENOTYPE_MUTATION = pheno_mutations[xx];

for (int yy = 0; yy < num_strats; yy++)
{

        STRATEGY_MUTATION = strat_mutations[yy];

for (int zz = 0; zz < num_benefits; zz++)
{

        BENEFIT = benefits[zz];

for (int xxx = 0; xxx < num_thresholds1; xxx++)
{

        HARDCODEDTHRESHOLD2 = larger_thresholds[xxx];


for (int yyyy = 0; yyyy < xxx+1; yyyy++)
{

	HARDCODEDTHRESHOLD1 = smaller_thresholds[yyyy];








/*Set the initial phenotype and strategy distribution */

	int time = 0;

        for(int i = 0; i < POPULATION; i++)
        {

                for(int j = 0; j < PHENOTYPES; j++)
                {

                        members[i][j] = 300;
                }

		members[i][PHENOTYPES] = INITIAL_POPULATION;
	
                results[i] = 0;

        }



/*This is just a bunch of places to store statistical meaures about the data to print later. */

int summing_coop = 0;
int total_coop = 0;

int summing_defec = 0;
int total_defec = 0;

int mutual_coop = 0;
int mutual_defec = 0;
int duped = 0;

double coop_variance = 0;
double defec_variance = 0;
double total_variance = 0;

double coop_mean = 0;
double defec_mean = 0;
double total_mean = 0;


int sum_dummy1 = 0;
int sum_dummy2 = 0;


/*Making a file to put all the data into */


int dummy1;	
char fnametemplate[256];
FILE *data_output;
snprintf(fnametemplate, 256, "Invasion_Pop%d_Gen%d_SS%f_PM%f_SM%f_IC%d_BURN%d_B%f_SMALLTHRES%f_SMALLTHRES%f_XXX", POPULATION, GENERATIONS, SELECTIONSTRENGTH, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, HARDCODEDTHRESHOLD1, HARDCODEDTHRESHOLD2);
dummy1 = mkstemp(fnametemplate);
data_output = fdopen(dummy1, "w");


/*Header to the data */

fprintf(data_output, "STR_MUT, PHE_MUT, B, THRES1, THRES2, time, COOPs, DEFs, COOP_VAR, DEF_VAR, TOTAL_VAR, COOP_MEAN, DEF_MEAN, TOTAL_MEAN, CC_GAMES, DD_GAMES, CD_GAMES\n");



/*Here is the actual time loop */



        for(int i = 0; i < GENERATIONS; i++)
        {

/*NOTICE: WHEN YOU ADD MORE TYPES OF PHENOTYPES YOU WILL NEED TO CHANGE THIS PRINTING, ITS JUST FOR ONE PHENOTYPE. */

/*We need to print the phenotypes to make a dancing phenotype histogram GIF. This prints directly to the pipeline, but the seperate file
prints the relevant statistical data. */


printf("%d \n", i);

		for (int j = 0; j < POPULATION; j++)
		{

			if(members[j][PHENOTYPES] == 1)
			{

				printf("%lf, ", members[j][PHENOTYPES-1]);
				total_coop += 1;

			}
		}
	
		printf("\n");

		for (int j = 0; j < POPULATION; j++)
		{

			if(members[j][PHENOTYPES] == -1)
			{

				printf("%lf, ", members[j][PHENOTYPES - 1]);
				total_defec += 1;

			}
		}

		
		printf("\n");

/*This is to get the means and variances and so forth of the cooperators/defectors. */

		double coop_pheno[total_coop];
		double defec_pheno[total_defec];
		double total_pheno[total_coop+total_defec];

		
		for (int j = 0; j < POPULATION; j++)
		{

			if(members[j][PHENOTYPES] == 1)
                        {

				coop_pheno[sum_dummy1] = members[j][PHENOTYPES - 1];
				sum_dummy1 += 1;

			}

			if(members[j][PHENOTYPES] == -1)
			{

				defec_pheno[sum_dummy2] = members[j][PHENOTYPES - 1];
				sum_dummy2 += 1;

			}

			total_pheno[j] = members[j][PHENOTYPES - 1];		

		}

		sum_dummy1 = 0;
		sum_dummy2 = 0;

		

/* If the file has -50 or N/A, there's an issue, such as being unable to take the variance of 1 value. */

double coop_variance = -50;
double defec_variance = -50;
double total_variance = -50;

double coop_mean = -50;
double defec_mean = -50;
double total_mean = -50;

		if(total_coop != 0)
		{
			coop_variance = gsl_stats_variance(coop_pheno, 1, total_coop);
			coop_mean = gsl_stats_mean(coop_pheno, 1, total_coop);
		}

		if(total_defec != 0)
		{
			defec_variance = gsl_stats_variance(defec_pheno, 1, total_defec);
			defec_mean = gsl_stats_mean(defec_pheno, 1, total_defec);
		}

		total_variance = gsl_stats_variance(total_pheno, 1, (total_defec+total_coop));

		total_mean = gsl_stats_mean(total_pheno, 1, (total_defec+total_coop));


		interaction(members, results, time, &mutual_coop, &mutual_defec, &duped);

                time = i;
		
	
                fprintf(data_output, "%f, %f, %f, %f, %f, %d, %d, %d, %f, %f, %f, %f, %f, %f, %d, %d, %d\n", STRATEGY_MUTATION, PHENOTYPE_MUTATION, BENEFIT, HARDCODEDTHRESHOLD1, HARDCODEDTHRESHOLD2, time, total_coop, total_defec, coop_variance, defec_variance, total_variance, coop_mean, defec_mean, total_mean, mutual_coop, mutual_defec, duped);

//printf("%f, %f, %f, %f, %f, %d, %d, %d, %f, %f, %f, %f, %f, %f, %d, %d, %d\n", STRATEGY_MUTATION, PHENOTYPE_MUTATION, BENEFIT, HARDCODEDTHRESHOLD1, HARDCODEDTHRESHOLD2, time, total_coop, total_defec, coop_variance, defec_variance, total_variance, coop_mean, defec_mean, total_mean, mutual_coop, mutual_defec, duped);

                total_coop = 0;
                total_defec = 0;

                newgeneration(r, members, results, time);

	}

	fclose(data_output);


}
}
}
}
}




}

void interaction(double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, int *mutual_coop, int *mutual_defec, int *duped)
{

*mutual_coop = 0;
*mutual_defec = 0;
*duped = 0;


double distance = 0;

/*This is all the games that occur. We go through each individual, and they play a game with every other individual. 

NOTE::: HARD CODED THRESHOLD 2 MUST BE SPECIFIED AS THE LARGER value.*/

	for (int i = 0; i < POPULATION; i++)
	{

		for (int j = 0; j < POPULATION; j++)
		{

			if(i != j)
			{
			
				distance = 0;

				for(int k = 0; k < PHENOTYPES; k++)
				{

					distance += (members[i][k] - members[j][k])*(members[i][k] - members[j][k]);
				}

				
				if (sqrt(distance) > HARDCODEDTHRESHOLD2)
				{

					results[i] += COST;
					*mutual_defec += 1;

				}

				if (sqrt(distance) < HARDCODEDTHRESHOLD1)
				{

					results[i] += BENEFIT;
					*mutual_coop += 1;

				}

				if (sqrt(distance) > HARDCODEDTHRESHOLD1 && sqrt(distance) < HARDCODEDTHRESHOLD2)
				{


					if (members[i][PHENOTYPES] == -1 && members[j][PHENOTYPES] == 1)
					{

						results[i] += BENEFIT + COST;
						*duped += 1;
					}

					else if (members[i][PHENOTYPES] == 1 && members[j][PHENOTYPES] == -1)
					{
						results[i] += 0;
						*duped += 1;
					}

					else if(members[i][PHENOTYPES] == 1 && members[j][PHENOTYPES] == 1)
					{
						results[i] += BENEFIT;
						*mutual_coop += 1;
					}
					else if(members[i][PHENOTYPES] == -1 && members[j][PHENOTYPES] == -1)
					{

						results[i] += COST;
						*mutual_defec += 1;

					}
				}


			}

		}

	}

}



void newgeneration(const gsl_rng * r, double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time)
{

/*Now we repopulate for the next generation, weighted by the success in the previous games.*/

double fitness_weights[POPULATION];

	for (int i = 0; i < POPULATION; i++)
	{

		fitness_weights[i] = 1 + (results[i])*SELECTIONSTRENGTH;

	}

	double shadow_members[POPULATION][PHENOTYPES+1];

	for (int i = 0; i < POPULATION; i++)
	{

		for (int j = 0; j < PHENOTYPES+1; j++)
		{

			shadow_members[i][j] = members[i][j];

		}

	}

/*Use the multinomial draw from GSL */

	int K = POPULATION; /* # "categories" */
	unsigned int n[K]; /*array to return results of draw */

	gsl_ran_multinomial(r, K, POPULATION, fitness_weights, n);
	
	int counter = 0;

/*Here we fill the array with all the individuals listed in n[] */

	for (int i = 0; i < POPULATION; i++)
	{

		for(int j = 0; j < n[i]; j++)
		{

			members[j+counter][0] = shadow_members[i][0];
			members[j+counter][PHENOTYPES] = shadow_members[i][PHENOTYPES];

		}

		counter += n[i];

	}

int random_dummy = 0;

/*add a phenotype mutation */

	for (int i = 0; i < POPULATION; i++)
	{

		for(int j = 0; j < PHENOTYPES; j++)
		{

			members[i][j] += gsl_ran_gaussian(r, PHENOTYPE_MUTATION);

		}

	}



/* If time is greater than the burn in, we allow strategy mutation. This is the "invasion" aspect.*/

	if (time > BURN_IN)
	{
		for (int i = 0; i < POPULATION; i++)
		{

			if(gsl_ran_bernoulli(r, STRATEGY_MUTATION))
			{

				
					members[i][PHENOTYPES] *= -1;
			}

		}
	}

	for (int i = 0 ; i < POPULATION; i++)
	{

		results[i] = 0;
	}


}





