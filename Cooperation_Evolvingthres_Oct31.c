/* This code was originally written by Linnea Bavik during her Spring 2019 research rotation at Emory University in the Weissman lab. It describes a model for the evolution of phenotype-based cooperation strategies in a population.*/

/*August 7 update: Any dimension of continuous phenotype is possible. Only two manually input strategies are possible.

The program can output:

(i) Statistics file. This is a csv file with multiple columns labelled as such:

STR_MUT, PHE_MUT, B, STRICT_THRES, 
GENEROUS_THRES, time, STRICTs, GENEROUSs, CC_GAMES, DD_GAMES, CD_GAMES, GENEROUS_MEAN_%d, GENEROUS_VAR_%d, STRICT_MEAN_%d, STRICT_VAR_%d, TOTAL_MEAN_%d, TOTAL_VAR_%d

Here %d refers to 1, ..., number of specified phenotype dimensions. So the header could list: 

..., 1_STRICT_MEAN, 1_GENEROUS_MEAN, 1_STRICT_VAR, 1_GENEROUS_VAR, 2_STRICT_MEAN, 2_GENEROUS_MEAN, 2_STRICT_VAR, 2_GENEROUS_VAR, ... etc.


(ii) Phenotypes file. Each row is a list of phenotypes of a given threshold. Best served for histogram plotting in Mathematica. Same information is given in the Fitnesses file (v).


(iii) Parent Map file. This gives the mutations and their current number of surviving descendants.


(iv) Radial Density file. This is the counts for a histogram of all the distances from the center of mass for each memebr of the population. Same information is given in the Fitnesses file (v).


(v) Fitnesses file. A csv file with columns labelled:

PHENOTYPE_%d, ..., TIME, INDIVIDUAL, DISTANCE, THRESHOLD, SCORE


(vi) Pairwise Distances file. Each row is a list of pairwise distances, first between Generous-Generous, then Strict-Strict, then Strict-Generous. This information can be derived from the output of the Fitnesses file, using data processing in Python or the like. (Typically this is what we've done to cut down on extraneous big files hanging around.)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>


/*
This game is a shifted Donation Game. The game matrix gives benefit B for mutual cooperation, benefit C for mutual defection, B+C for defecting upon one's partner and 0 for being defected on.
*/

int ITERATIONS = 1;

int POPULATION = 1000;
int GENERATIONS = 5000;
int PHENOTYPES = 10;
double SELECTIONSTRENGTH = 1;
double PHENOTYPE_MUTATION = 1;
double STRATEGY_MUTATION = 0.0001;

double ZETA = -0.001;

int INITIAL_POPULATION = 3;
int BURN_IN = 0;

double COST = 1;
double BENEFIT = 10;

double GENEROUS_THRESHOLD = 10;
double STRICT_THRESHOLD = 5;



void interaction(double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, int *mutual_coop, int *mutual_defec, int *duped);

void newgeneration(const gsl_rng * r, double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, double strategy_map[POPULATION], double shadow_strategy_map[POPULATION], double strategy_origins[1000], int *strategy_increment);

#define num_phenos 1
#define num_strats 4
#define num_benefits 1
#define num_generous_thresholds 1
#define num_strict_thresholds 1
#define num_zetas 2

/* These are arrays of values one wishes to test. This is very easily parallelizable, but hasn't been done yet. */


double pheno_mutations[num_phenos] = {1};
double strat_mutations[num_strats] = {0.001, 0.01, 0.1, 1};
double benefits[num_benefits] = {10};
double generous_thresholds[num_generous_thresholds] = {10};
double strict_thresholds[num_strict_thresholds] = {1};
double zetas[num_zetas] = {-0.001, -0.01};


/*the following variable declared for the parameter sweep */

#define num_precision_ranges 1

int main()
{

/* Here we set which outputs are desired, set true for desired files, set step_size for skips between outputs (i.e. do you want every time step (step_size = 1) ? Or every five time steps? */

	bool stat_on = true;
	bool fitnesses_on = true;
	bool pheno_on = true;
	bool radial_on = false;
	bool parent_on = false;
	bool pairs_on = false;

	int step_size = 1;


/* For adjusting parameter sweep values */

/*For GSL Library multinomial draw */

	const gsl_rng_type *T; 
/*generator type */

	gsl_rng * r; 
/*creates a random number generator */

	gsl_rng_env_setup(); 
/*read from environment variable */

	gsl_rng_default_seed = time(NULL);

/* choose defult generator type which I set to the mersenne twister. */
	gsl_rng_default = gsl_rng_mt19937;

	T = gsl_rng_default;

/*creates an instance */
	r = gsl_rng_alloc (T); 

        srand(time(NULL));



        double members[POPULATION][PHENOTYPES+1];
        double results[POPULATION];


		double strategy_map[POPULATION];
		double shadow_strategy_map[POPULATION];
		double strategy_origins[1200];
		int strategy_increment = 1;


/*HERE's THE BIG OL LOOP */

for (int pheno_x = 0; pheno_x < num_phenos; pheno_x++)
{

        PHENOTYPE_MUTATION = pheno_mutations[pheno_x];

for (int stratmut_x = 0; stratmut_x < num_strats; stratmut_x++)
{

        STRATEGY_MUTATION = strat_mutations[stratmut_x];

for (int benefit_x = 0; benefit_x < num_benefits; benefit_x++)
{

        BENEFIT = benefits[benefit_x];

for (int zeta_x = 0; zeta_x < num_zetas; zeta_x++)
{

	ZETA = zetas[zeta_x];

for (int iter = 1; iter < ITERATIONS+1; iter++)
{


/*Set the initial phenotype and strategy distribution */

	int time = 0;

        for(int i = 0; i < POPULATION; i++)
        {

                for(int j = 0; j < PHENOTYPES; j++)
                {

                        members[i][j] = 300;
                }

		members[i][PHENOTYPES] = INITIAL_POPULATION;

		if(INITIAL_POPULATION == 3)
		{

			members[i][PHENOTYPES] = rand()%50;
//			members[i][PHENOTYPES] += gsl_ran_gaussian(r, PHENOTYPE_MUTATION);

		}

		if(parent_on == true)
		{
			strategy_map[i] = 1;
		}
	
                results[i] = 0;

        }



/*This is just a bunch of places to store statistical meaures about the data to print later. */

	int mutual_coop = 0;
	int mutual_defec = 0;
	int duped = 0;

	double total_variance = 0;

	double total_skew = 0;

	double total_kurtosis = 0;
	
	double total_mean = 0;

	double total_pheno[POPULATION];

	


	int sum_dummy1 = 0;
	int sum_dummy2 = 0;

	FILE *stat_output, *fitness_output, *pheno_output, *radial_output, *parent_output, *pairs_output;

/*Creating the files to output to. */

	if(stat_on == true)
	{

		int dummy_stat;
		char stat_template[256];
//		FILE *stat_output;

		snprintf(stat_template, 256, "STATISTICS_EVOLTHRES_Dim%.2d_Pop%.2d_Gen%.2d_PM%.5f_SM%.5f_IC%.2d_BURN%.2d_B%.2f_ZETA%.5f_iteration%d", PHENOTYPES, POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, ZETA, iter);

		stat_output = fopen(stat_template, "w");

	}

	if(fitnesses_on == true)
        {

                int dummy_fitness;
                char fitness_template[256];
//                FILE *fitness_output;

                snprintf(fitness_template, 256, "FITNESSES_EVOLTHRES_Dim%.2d_Pop%.2d_Gen%.2d_PM%.5f_SM%.5f_IC%.2d_BURN%.2d_B%.2f_ZETA%.5f_iteration%d", PHENOTYPES, POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, ZETA, iter);


//                dummy_fitness = mkstemp(fitness_template);

//              fitness_output = fdopen(dummy_fitness, "w");

		fitness_output = fopen(fitness_template, "w");

        }

	if(pheno_on == true)
        {

                int dummy_pheno;
                char pheno_template[256];
//              FILE *pheno_output;

                snprintf(pheno_template, 256, "PHENOTYPES_EVOLTHRES_Dim%.2d_Pop%.2d_Gen%.2d_PM%.5f_SM%.5f_IC%.2d_BURN%.2d_B%.2f_ZETA%.5f_iteration%d", PHENOTYPES, POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, ZETA, iter);


//                dummy_pheno = mkstemp(pheno_template);

                pheno_output = fopen(pheno_template, "w");

        }

	if(radial_on == true)
        {

                int dummy_radial;
                char radial_template[256];
//                FILE *radial_output;

                snprintf(radial_template, 256, "RADIALDENSITY_EVOLTHRES_Dim%.2d_Pop%.2d_Gen%.2d_PM%.5f_SM%.5f_IC%.2d_BURN%.2d_B%.2f_ZETA%.5f_iteration%d", PHENOTYPES, POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, ZETA, iter);


                radial_output = fopen(radial_template, "w");

        }

	
	if(pairs_on == true)
        {

                int dummy_pairs;
                char pairs_template[256];
//                FILE *pairs_output;

                snprintf(pairs_template, 256, "PAIRWISE_EVOLTHRES_Dim%.2d_Pop%.2d_Gen%.2d_PM%.5f_SM%.5f_IC%.2d_BURN%.2d_B%.2f_ZETA%.5f_iteration%d", PHENOTYPES, POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, ZETA, iter);


                pairs_output = fopen(pairs_template, "w");

        }


/*Header to the stat data */

	if(stat_on == true)
	{
		fprintf(stat_output, "STR_MUT, PHE_MUT, B, time, CC_GAMES, DD_GAMES, CD_GAMES");

	

		for (int i = 0; i < PHENOTYPES; i++)
		{

			if(stat_on == true)
			{
				fprintf(stat_output, ", PHENO_MEAN_%d, PHENO_VAR_%d, PHENO_SKEW_%d, PHENO_KURT_%d", (i+1), (i+1), (i+1), (i+1));
			}

			if(fitnesses_on == true)
			{
				fprintf(fitness_output, "PHENOTYPE_%d, ", i+1);
			}

		}

		fprintf(stat_output, ", THRES_MEAN, THRES_VAR, THRES_SKEW, THRES_KURTOSIS");

	}

	if(fitnesses_on == true)
	{
		fprintf(fitness_output, "TIME, INDIVIDUAL, DISTANCE, THRESHOLD, SCORE\n");
	}

	if(stat_on == true)
	{

		fprintf(stat_output, "\n");
	}


	double radius[POPULATION][POPULATION];

	if(radial_on == true)
	{

		for (int i = 0; i < POPULATION; i++)
		{

			for (int j = 0; j < POPULATION; j++)
			{

				radius[i][j] = 0;

			}

		}
	}


	double total_mean_vector[PHENOTYPES];

	double distance = 0;
	
	double threshold = 0;

	for (int inty = 0; inty < PHENOTYPES; inty++)
	{
		total_mean_vector[inty] = 0;
	}


/*Here is the actual time loop */

	int first = 0;

        for(int i = 0; i < GENERATIONS; i++)
        {
		time = i;

		interaction(members, results, time, &mutual_coop, &mutual_defec, &duped);

		if (fitnesses_on == true && i%step_size == 0)
		{

			for (int iterator = 0; iterator < POPULATION; iterator++)
			{

				for (int free_int = 0; free_int < PHENOTYPES; free_int++)
				{

					distance += (members[iterator][free_int] - total_mean_vector[free_int])*(members[iterator][free_int] - total_mean_vector[free_int]);

					fprintf(fitness_output, "%.3f, ", members[iterator][free_int]);

				}			

				distance = sqrt(distance);

				threshold = members[iterator][PHENOTYPES];

				fprintf(fitness_output, "%d, %d, %.3f, %.3f, %.5f\n", i, iterator, distance, threshold, results[iterator]);

			}


		}


                newgeneration(r, members, results, time, strategy_map, shadow_strategy_map, strategy_origins, &strategy_increment);


/*Output radial distances: */
		

		if(radial_on == true && i%step_size == 0)
		{
		
			for (int j = 0; j < PHENOTYPES; j++)
			{

				for (int k = 0; k < POPULATION; k++)
				{

					for (int l = k; l < POPULATION; l++)
					{

						radius[k][l] += fabs(members[k][j] - members[l][j])*fabs(members[k][j] - members[l][j]);

					}

				}

			}


			for (int j = 0; j < POPULATION; j++)
			{

				for (int k = j; k < POPULATION; k++)
				{

					if(first == 0)
					{

						fprintf(radial_output, "%lf", sqrt(radius[j][k]));
						first = 1;
					}	

					else
					{

						fprintf(radial_output, ", %lf", sqrt(radius[j][k]));

					}

				}

			}

			first = 0;

			fprintf(radial_output, "\n");
		}

/*Output strategies from the new generation for the parent map*/

		if(parent_on == true && i%step_size == 0)
		{
                
			for (int j = 0; j < POPULATION; j++)
                	{

				if (first == 0)
				{

					fprintf(parent_output, "%f", strategy_map[j]);
					first = 1;
				}
			
				else
				{
	                        	fprintf(parent_output, ", %f", strategy_map[j]);

				}

	                }

			first = 0;

			fprintf(parent_output, "\n");
		}


/* Output phenotypes for the phenotype distribution */

		if(pheno_on == true && i%step_size == 0)
		{


			for (int j = 0; j < PHENOTYPES + 1; j++)
			{


				if(pheno_on == true && i%step_size == 0)
				{
					fprintf(pheno_output, "PHENOTYPE %d: (if number > pheno dim, its the thresholds)\n", (j+1));
				}

				for (int k = 0; k < POPULATION; k++)
				{

/* SPACE SAVE CONDITIONS */
					if(pheno_on == true && i%step_size == 0)
					{
					
						if (first == 0)
						{
	
							fprintf(pheno_output, "%lf", members[k][j]);

							first = 1;

						}

						else
						{

							fprintf(pheno_output, ", %lf", members[k][j]);
						}
					}
/* I set j == 0 here because we only want to count each loose individual once */

				}

			}

			first = 0;

			if(pheno_on == true && i%step_size == 0)
			{

				fprintf(pheno_output, "\n");			

			}
			
		}

/* Here we add the pairwise stuff: */

		if (pairs_on == true && i%step_size == 0)
		{

			double p_dist = 0.0;

			for (int p_i = 0; p_i < POPULATION; p_i++)
			{

				for (int p_j = p_i + 1; p_j < POPULATION; p_j++)
				{

					for(int p_k = 0; p_k < PHENOTYPES; p_k++)
					{

						p_dist += (members[p_i][p_k] - members[p_j][p_k])*(members[p_i][p_k] - members[p_j][p_k]);

					}

					fprintf(pairs_output, "%.3f, ", sqrt(p_dist));

					p_dist = 0.0;

				}			

			}

			fprintf(pairs_output, "\n");

		}


		if(stat_on == true && i%step_size == 0)
		{
			fprintf(stat_output, "%f, %f, %f, %d, %d, %d, %d", STRATEGY_MUTATION, PHENOTYPE_MUTATION, BENEFIT, time, mutual_coop, mutual_defec, duped);
		}

/*Now to get the means and variance of the phenotypes. */

		total_variance = -50;
		total_mean = -50;
		total_skew = -50;
		total_kurtosis = -50;

		sum_dummy1 = 0;
		sum_dummy2 = 0;

		for (int j = 0; j < PHENOTYPES+1; j++)
		{

			for (int k = 0; k < POPULATION; k++)
			{

				total_pheno[k] = members[k][j];

			}

			sum_dummy1 = 0;
			sum_dummy2 = 0;

			total_variance = gsl_stats_variance(total_pheno, 1, POPULATION);

			total_mean = gsl_stats_mean(total_pheno, 1, POPULATION);

			total_skew = gsl_stats_skew(total_pheno, 1, POPULATION);

			total_kurtosis = gsl_stats_kurtosis(total_pheno, 1, POPULATION);

			total_mean_vector[j] = total_mean;

			if(stat_on == true && i%step_size == 0)
			{	
				fprintf(stat_output, ", %f, %f, %f, %f", total_mean, total_variance, total_skew, total_kurtosis);
			}

		}

		if (stat_on == true && i%step_size == 0)
		{
			fprintf(stat_output, "\n");
		}
/*
		for (int countt = 0; countt < POPULATION; countt++)
		{
			printf("%f, ", members[countt][PHENOTYPES]);
		}

		printf("\n");
*/
	}

	if(radial_on == true)
	{
		fclose(radial_output);
	}

	if(pheno_on == true)
	{
		fclose(pheno_output);
	}

	if(stat_on == true)
	{
		fclose(stat_output);
	}

	if(fitnesses_on == true)
	{
		fclose(fitness_output);
	}

	int i = 0;

	if(parent_on == true)
	{
		fprintf(parent_output, "origins:\n");
/*
		while(strategy_origins[i] < GENERATIONS - 5)
		{

			fprintf(parent_output, "%f, ", strategy_origins[i]);
			i++;

		} 
*/

		fclose(parent_output);

	}

printf("\n \n \n \n");

}
}

}
}
}
}
//}
//}













void interaction(double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, int *mutual_coop, int *mutual_defec, int *duped)
{

*mutual_coop = 0;
*mutual_defec = 0;
*duped = 0;


double distance = 0;

/*This is all the games that occur. We go through each individual, and they play a game with every other individual. 
NOTE::: HARD CODED THRESHOLD 2 MUST BE SPECIFIED AS THE LARGER value.*/

	for (int i = 0 ; i < POPULATION; i++)
        {

                results[i] = 0;
        }



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

/*EVOLVING THRESHOLDS BIT */

				if(members[i][PHENOTYPES] == members[j][PHENOTYPES])
				{

					if(distance > members[i][PHENOTYPES])
					{

						results[i] += COST;
						*mutual_defec += 1;

					} else
					{

						results[i] += BENEFIT;
						*mutual_coop += 1;

					}

				}

				else
				{


					GENEROUS_THRESHOLD = fmax(members[i][PHENOTYPES], members[j][PHENOTYPES]);
					STRICT_THRESHOLD = fmin(members[i][PHENOTYPES], members[j][PHENOTYPES]);
				
				
					if (sqrt(distance) > GENEROUS_THRESHOLD)
					{

						results[i] += COST;
						*mutual_defec += 1;

					}

					if (sqrt(distance) < STRICT_THRESHOLD)
					{

						results[i] += BENEFIT;
						*mutual_coop += 1;

					}

					if (sqrt(distance) > STRICT_THRESHOLD && sqrt(distance) < GENEROUS_THRESHOLD)
					{
	

						if (members[i][PHENOTYPES] == STRICT_THRESHOLD)
						{

							results[i] += BENEFIT + COST;
							*duped += 1;
						}

						else if (members[i][PHENOTYPES] == GENEROUS_THRESHOLD)
						{
							results[i] += 0;
							*duped += 1;
						}

					}

				}

			}

		}

	}

}





void newgeneration(const gsl_rng * r, double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, double strategy_map[POPULATION], double shadow_strategy_map[POPULATION], double strategy_origins[1000], int *strategy_increment)
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
/*PARENT MAP STUFF ~~~~~~~~~*/


		shadow_strategy_map[i] = strategy_map[i];

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

			for (int k = 0; k <= PHENOTYPES; k++)
			{
				members[j+counter][k] = shadow_members[i][k];
			}

/*PARENT MAP STUFF ~~~~~~~~~~~~~~~~~ */
			strategy_map[j+counter] = shadow_strategy_map[i];

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

//			members[i][PHENOTYPES] *= gsl_ran_lognormal(r, ZETA, sqrt(-1*2*ZETA));

			if(gsl_ran_bernoulli(r, STRATEGY_MUTATION))
			{

				members[i][PHENOTYPES] *= gsl_ran_lognormal(r, ZETA, sqrt(-1*2*ZETA));
				
//					members[i][PHENOTYPES] *= -1;
/*PARENT MAP STUFF*/	
//printf("mapping %f to %f \n", strategy_map[i], (-1)*strategy_map[i]/(fabs(strategy_map[i]))*(fabs(strategy_map[i]) + 1));

					strategy_map[i] = (-1)*strategy_map[i]/(fabs(strategy_map[i]))*(*strategy_increment + 1);

//printf("strategy_origins at *strategy increment (%d) is %f before \n", *strategy_increment, strategy_origins[*strategy_increment]);

//					strategy_origins[*strategy_increment] = time;

//printf("strategy_origins at *strategy increment (%d) is %f after\n", *strategy_increment, strategy_origins[*strategy_increment]);

					*strategy_increment += 1;


			}

		}
	}

/*	for (int i = 0 ; i < POPULATION; i++)
	{

		results[i] = 0;
	}
*/

}



