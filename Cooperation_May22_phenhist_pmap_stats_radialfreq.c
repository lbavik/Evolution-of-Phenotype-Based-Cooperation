/* This code was written by Linnea Bavik during her Spring 2019 research rotation at Emory University in the Weissman lab.
It describes a model for the evolution of phenotype-based cooperation strategies in a population.*/

/* May 22 update: Multiple phenotypes are possible. The program outputs four outputs: (i) phenotypes of every member of the population at every time step for histogram plotting. (ii) a parent map: gives the mutations and their current number of surviving descendants. (iii) the statistics which is a list under the header:

STR_MUT, PHE_MUT, B, STRICT_THRES, LOOSE_THRES, time, STRICTs, LOOSEs, CC_GAMES, DD_GAMES, CD_GAMES, LOOSE_MEAN_%d, LOOSE_VAR_%d, STRICT_MEAN_%d, STRICT_VAR_%d, TOTAL_MEAN_%d, TOTAL_VAR_%d

Here %d refers to 1, ..., number of specified phenotype dimensions. So the header could list: 

..., 1_STRICT_MEAN, 1_LOOSE_MEAN, 1_STRICT_VAR, 1_LOOSE_VAR, 2_STRICT_MEAN, 2_LOOSE_MEAN, 2_STRICT_VAR, 2_LOOSE_VAR, ... etc.

(iv) radial density. (counts for a histogram of all the distances from each member of the population in the multidimensional phenotype space.)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>


/*
The May 22 2019 version of this code refers to population evolution in an invasion scenario when individuals use either a strict or
loose hard threshold of phenotype distance to decide cooperation or defection. The game matrix gives benefit B for mutual cooperation,
benefit C for mutual defection, B+C for defecting upon one's partner and 0 for being defected on.
*/



int POPULATION = 1000;
int GENERATIONS = 5000;
int PHENOTYPES = 10;
double SELECTIONSTRENGTH = 1;
double PHENOTYPE_MUTATION = 1;
double STRATEGY_MUTATION = .02;

int INITIAL_POPULATION = -1;
int BURN_IN = 70000;

double COST = 1;
double BENEFIT = 1.1;

double LOOSE_THRESHOLD = 10;
double STRICT_THRESHOLD = 5;



void interaction(double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, int *mutual_coop, int *mutual_defec, int *duped);

void newgeneration(const gsl_rng * r, double members[POPULATION][PHENOTYPES+1], double results[POPULATION], int time, double strategy_map[POPULATION], double shadow_strategy_map[POPULATION], double strategy_origins[1000], int *strategy_increment);

#define num_phenos 1
#define num_strats 1
#define num_benefits 3
#define num_loose_thresholds 1
#define num_strict_thresholds 23

/* These are arrays of values one wishes to test. This is very easily parallelizable, but hasn't been done yet. */


double pheno_mutations[num_phenos] = {1};
double strat_mutations[num_strats] = {.01};
double benefits[num_benefits] = {1, 2, 5};
double loose_thresholds[num_loose_thresholds] = {.1};
double strict_thresholds[num_strict_thresholds] = {0, .25, .5, .75, 1, 1.25, 1.50, 1.75, 2.0, 2.25, 2.50, 2.75, 3, 6, 9, 25, 50, 75, 100, 250, 500, 750, 1000};


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
/*PARENT MAP STUFF ~~~~~~~~~~~~~~~~~~~~~*/
	double strategy_map[POPULATION];
	double shadow_strategy_map[POPULATION];
	double strategy_origins[1200];
//	double strategy_origins[(int)(POPULATION*GENERATIONS*STRATEGY_MUTATION*1.1)];
	int strategy_increment = 1;



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

for (int xxx = 0; xxx < num_loose_thresholds; xxx++)
{

        LOOSE_THRESHOLD = loose_thresholds[xxx];


for (int yyyy = 0; yyyy < num_strict_thresholds; yyyy++)
{

	STRICT_THRESHOLD = strict_thresholds[yyyy];








/*Set the initial phenotype and strategy distribution */

	int time = 0;

        for(int i = 0; i < POPULATION; i++)
        {

                for(int j = 0; j < PHENOTYPES; j++)
                {

                        members[i][j] = 300;
                }

		members[i][PHENOTYPES] = INITIAL_POPULATION;

/*P MAP ~~~~~~~~~~~~~*/
		strategy_map[i] = 1;
		
	
                results[i] = 0;

        }



/*This is just a bunch of places to store statistical meaures about the data to print later. */

int summing_loose = 0;
int total_loose = 0;

int summing_strict = 0;
int total_strict = 0;

int mutual_coop = 0;
int mutual_defec = 0;
int duped = 0;

double loose_variance = 0;
double strict_variance = 0;
double total_variance = 0;

double loose_mean = 0;
double strict_mean = 0;
double total_mean = 0;


int sum_dummy1 = 0;
int sum_dummy2 = 0;


/*Making the four files to put all the data into: */


int dummy_pheno, dummy_parent, dummy_stat, dummy_radial;
	
char pheno_template[256], parent_template[256], stat_template[256], radial_template[256];
FILE *pheno_output, *parent_output, *stat_output, *radial_output;


snprintf(pheno_template, 256, "PHENOTYPES_Pop%.2d_Gen%.2d_PM%.2f_SM%.2f_IC%.2d_BURN%.2d_B%.2f_STRTHRES%.2f_LOOSETHRES%.2f", POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, STRICT_THRESHOLD, LOOSE_THRESHOLD);

snprintf(parent_template, 256, "PARENTMAP_Pop%.2d_Gen%.2d_PM%.2f_SM%.2f_IC%.2d_BURN%.2d_B%.2f_STRTHRES%.2f_LOOSETHRES%.2f", POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, STRICT_THRESHOLD, LOOSE_THRESHOLD);

snprintf(stat_template, 256, "STATISTICS_Pop%.2d_Gen%.2d_PM%.2f_SM%.2f_IC%.2d_BURN%.2d_B%.2f_STRTHRES%.2f_LOOSETHRES%.2f", POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, STRICT_THRESHOLD, LOOSE_THRESHOLD);

snprintf(radial_template, 256, "RADIALDENSITY_Pop%.2d_Gen%.2d_PM%.2f_SM%.2f_IC%.2d_BURN%.2d_B%.2f_STRTHRES%.2f_LOOSETHRES%.2f", POPULATION, GENERATIONS, PHENOTYPE_MUTATION, STRATEGY_MUTATION, INITIAL_POPULATION, BURN_IN, BENEFIT, STRICT_THRESHOLD, LOOSE_THRESHOLD);


dummy_pheno = mkstemp(pheno_template);
dummy_parent = mkstemp(parent_template);
dummy_stat = mkstemp(stat_template);
dummy_radial = mkstemp(radial_template);

pheno_output = fdopen(dummy_pheno, "w");
parent_output = fdopen(dummy_parent, "w");
stat_output = fdopen(dummy_stat, "w");
radial_output = fdopen(dummy_radial, "w");



/*Header to the stat data */

fprintf(stat_output, "STR_MUT, PHE_MUT, B, LOOSE_THRES, STRICT_THRES, time, LOOSEs, STRICTs, CC_GAMES, DD_GAMES, CD_GAMES");

	for (int i = 0; i < PHENOTYPES; i++)
	{

		fprintf(stat_output, ", LOOSE_MEAN_%d, LOOSE_VAR_%d, STRICT_MEAN_%d, STRICT_VAR_%d, TOTAL_MEAN_%d, TOTAL_VAR_%d", (i+1), (i+1), (i+1), (i+1), (i+1), (i+1));

	}

	fprintf(stat_output, "\n");



	double radius[POPULATION][POPULATION];

	for (int i = 0; i < POPULATION; i++)
	{

		for (int j = 0; j < POPULATION; j++)
		{

			radius[i][j] = 0;

		}

	}



/*Here is the actual time loop */

	int first = 0;

        for(int i = 0; i < GENERATIONS; i++)
        {

		time = i;

		interaction(members, results, time, &mutual_coop, &mutual_defec, &duped);

		newgeneration(r, members, results, time, strategy_map, shadow_strategy_map, strategy_origins, &strategy_increment);

/*Output radial distances: */
/*CONDITIONS TO SAVE SPACE~~~*/
if(i > GENERATIONS-50)
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
/*SPACE SAVE CONDITIONS */
/*
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
*/


/* Output phenotypes for the phenotype distribution */


		for (int j = 0; j < PHENOTYPES; j++)
		{
/*SPACE SAVE*/
if(i > GENERATIONS - 50)
{
			fprintf(pheno_output, "PHENOTYPE %d: (loose then strict)\n", (j+1));
}
			for (int k = 0; k < POPULATION; k++)
			{

				if(members[k][PHENOTYPES] == 1)
				{
/* SPACE SAVE CONDITIONS */
if(i > GENERATIONS-50)
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

					if(j == 0)
					{

						total_loose += 1;

					}


				}


			}

			first = 0;


			fprintf(pheno_output, "\n");			

			
			for (int k = 0; k < POPULATION; k++)
			{

				if(members[k][PHENOTYPES] == -1)
                                {
if(i > GENERATIONS - 50)
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
/* I set j == 0 here because we only want to count each strict individual once */

                                        if(j == 0)
                                        {

                                                total_strict += 1;

                                        }


                                }


			}

			first = 0;

			fprintf(pheno_output, "\n");

		}




		fprintf(stat_output, "%f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %d", STRATEGY_MUTATION, PHENOTYPE_MUTATION, BENEFIT, LOOSE_THRESHOLD, STRICT_THRESHOLD, time, total_loose, total_strict, mutual_coop, mutual_defec, duped);

/*Now to get the means and variance of the phenotypes. */

		loose_variance = -50;
		strict_variance = -50;
		total_variance = -50;

		loose_mean = -50;
		strict_mean = -50;
		total_mean = -50;

		double loose_pheno[total_loose];
		double strict_pheno[total_strict];
		double total_pheno[total_loose+total_strict];

		sum_dummy1 = 0;
		sum_dummy2 = 0;

		for (int j = 0; j < PHENOTYPES; j++)
		{

			for (int k = 0; k < POPULATION; k++)
			{

				if(members[k][PHENOTYPES] == 1)
				{

					loose_pheno[sum_dummy1] = members[k][j];
					sum_dummy1 += 1;

				}

				if(members[k][PHENOTYPES] == -1)
				{

					strict_pheno[sum_dummy2] = members[k][j];
					sum_dummy2 += 1;

				}


				total_pheno[k] = members[k][j];

			}

			sum_dummy1 = 0;
			sum_dummy2 = 0;

			if (total_loose != 0)
			{

				loose_variance = gsl_stats_variance(loose_pheno, 1, total_loose);

				loose_mean = gsl_stats_mean(loose_pheno, 1, total_loose);

			}


			if (total_strict != 0)
			{

				strict_variance = gsl_stats_variance(strict_pheno, 1, total_strict);
                                strict_mean = gsl_stats_mean(strict_pheno, 1, total_strict);

			}

			total_variance = gsl_stats_variance(total_pheno, 1, (total_strict+total_loose));

			total_mean = gsl_stats_mean(total_pheno, 1, (total_strict+total_loose));

			
			fprintf(stat_output, ", %f, %f, %f, %f, %f, %f", loose_mean, loose_variance, strict_mean, strict_variance, total_mean, total_variance);


		}

		total_loose = 0;
		total_strict = 0;

		fprintf(stat_output, "\n");

	}

	fclose(radial_output);
	fclose(pheno_output);
	fclose(stat_output);

	int i = 0;
/*
	fprintf(parent_output, "origins:\n");

	while(strategy_origins[i] < GENERATIONS - 5)
	{

		fprintf(parent_output, "%f, ", strategy_origins[i]);
		i++;

	} 
*/

	fclose(parent_output);

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

				
				if (sqrt(distance) > LOOSE_THRESHOLD)
				{

					results[i] += COST;
					*mutual_defec += 1;

				}

				if (sqrt(distance) < STRICT_THRESHOLD)
				{

					results[i] += BENEFIT;
					*mutual_coop += 1;

				}

				if (sqrt(distance) > STRICT_THRESHOLD && sqrt(distance) < LOOSE_THRESHOLD)
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

			if(gsl_ran_bernoulli(r, STRATEGY_MUTATION))
			{

				
					members[i][PHENOTYPES] *= -1;
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

	for (int i = 0 ; i < POPULATION; i++)
	{

		results[i] = 0;
	}


}



