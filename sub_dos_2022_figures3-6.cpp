/*
Created by Amanda Erin Wilson
Published to Github 08/26/2022

g++ -I eigen/ sub_dos_2022_figures3-6.cpp -o program2

Calculates expected distribution of states for subfunctionalization-only model and subfunctionalization + Dosage model after whole-genome duplication

To use this code need Eigen 
And to input:
// Nucleotide mutation rate, prob of mutations in coding region and reg region (lengths)
    long double nucleotide_mutation_rate;
    long double nucleotide_length_coding_region; //size of area where if a mutation occurs, that knocks out the coding region
    long double nucleotide_length_regulatory_region; //size of area where if a mutation occurs, that knocks out a regulatory region
// Number of regulatory regions and size of q matrix  
    int number_of_regulatory_regions;
    const int matrix_size; //this is equal to  number_of_regulatory_regions +2, is the size of the q matrix
    long double number_of_regulatory_regions_double;//this is number_of_regulatory_regions but as a double for calculation
// Thermodynamic variables (k, concentration of total A, concentration of total B)    
    long double k;
    long double a_total ;
    long double b_total ;     
    //[hp] = [A]free + [B]free
// Scalar on number of hydrophobic patches per cell and fitness penalty 
    long double scalar ;

// Figure parameters    
    // For naming .csv file
    std::string figure_number ;
    std::string trial_number ;    
    
    // Do you want to calculate for both Probability of fixation = Sella and Hirsch equation and 1/N?
    // Choose 1 if just for Sella and Hirsch, and 2 for both
    int number_of_prob_fix_calculations ;

    
    // effective Population sizes to loop through (N or Ne)
    std::vector <long double> pop_sizes ;
    int number_of_pops ;
    
    // How many time intervals to loop through (time at which we will calc distribution)
    int number_time_intervals ;
    // How big those time intervals should be
    long double size_time_intervals ; 
    

//AND change the initial state vector if need be

Eigen user manual: http://eigen.tuxfamily.org/dox/index.html


Pseudo code
    For each effective population size listed
        For each generation until desired time point
            Concentration of hydrophobic residues when in stoichiometric balance ([A]tot=[B]tot) = [A]free + [B]free (calculated from k, [A]total, [B]total)
            Concentration of hydrophobic residues when out of stoichiometric balance ([A]tot=2[B]tot) = [A]free + [B]free (calculated from k, [A]total, [B]total)

            For each state
                Unaffected tissues = z – current state
                Affected tissues  = current state
                Summation of hydrophobic residues across tissues = (unaffected tissues)*concentration hp when in stoichiometric balance + affected tissues * concentration hp when out of stoichiometric balance
                Fitness per state = 1/(scalar + sum of hp across tissues for state)
            Vector of fitness values for each state  = [fitness of state 0, fitness of state 1, etc]
            Probability Distribution Across States = using generator Markov matrix with vector of fitness values for each state, and initial state vector
*/
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include <ctime>
#include <string> 
#include <vector>
#include <fstream>


//------Function declaration-----------------------------------------   
std::string choose_which_prob_fix_string(int prob_fix_calculation);
long double calculate_concentration_in_bound_form(long double keq, long double x_total, long double y_total);
long double calculate_concentration_x_free(long double x_total, long double x_bound);
long double calculate_concentration_hydrophobic_residues(long double x_free, long double y_free);
long double do_all_stoichiometric_calc_for_concentration_hp(long double keq, long double x_total, long double y_total);
long double sum_hp_across_tissues(long double regulatory_regions, long double current_state, long double hp_per_tissue_in_balance, long double hp_per_tissue_out_of_balance);
long double calculate_fitness_of_state(long double sum_of_hp_across_tissues_for_state);
long double which_prob_fix_calculation(long double fitness_current_state, long double fitness_next_state, long double N, int prob_of_fix_which);
long double calculate_prob_fix_using_sella_hirsch(long double fitness_current_state, long double fitness_next_state, long double population_size);
long double calculate_prob_fix_as_one_over_Ne(long double population_size);
Eigen::MatrixXf calculate_probability_distribution_across_states(std::vector <long double> fitness_of_states, long double time, long double N, int prob_fix_calculation);

//------Initialization-----------------------------------------------

// Nucleotide mutation rate, prob of mutations in coding region and reg region (lengths)
    long double nucleotide_mutation_rate = 0.000000025;
    long double nucleotide_length_coding_region = 50000.00;
    long double nucleotide_length_regulatory_region = 775;
// Number of regulatory regions and size of q matrix  
    int number_of_regulatory_regions = 4;
    const int matrix_size = 6; //this is equal to  number_of_regulatory_regions +2
    long double number_of_regulatory_regions_double = 4.0;//this is number_of_regulatory_regions but as a double for calculation
// Thermodynamic variables (k, concentration of total A, concentration of total B)    
    long double k = 10000000000;
    long double a_total = 0.0000025;
    long double b_total = 0.0000025;     
// Scalar on number of hydrophobic patches per cell and fitness penalty 
    long double scalar = 1.00;

// Figure parameters    
    // For naming .csv file
    std::string figure_number = "3";
    std::string trial_number = "5000";    
    
    // Do you want to calculate for both Probability of fixation = Sella and Hirsch equation and 1/N?
    // Choose 1 if just for Sella and Hirsch, and 2 for both
    int number_of_prob_fix_calculations = 2;

    
    // Population sizes to loop through
    std::vector <long double> pop_sizes = {140000};
    int number_of_pops = 1;
    
    // How many time intervals to loop through
    int number_time_intervals = 1000;
    // How big those time intervals should be
    long double size_time_intervals = 5; 

//-------------------------------------------------------------------      
int main() {
// Name and Initiate .csv file to print to    
    std::ofstream myfile;
    std::string file_name = "fig" + figure_number + "_combined_vers4_wg" + trial_number +".csv";
    myfile.open (file_name);
// Label column headings    
    myfile << "Prob Fix,    Scalar, Ne, Keq,    Time,   ";
    for(int each_numbered_state = 0; each_numbered_state < number_of_regulatory_regions; each_numbered_state++){
        myfile << "State "<< each_numbered_state << ",   ";                    
    }
    myfile << "Subfuntionalized,    Pseudogenized \n";

// Loop through prob of fix = 1/N and Sella and Hirsch equation       
    for(int each_prob_fix_calculation = 0; each_prob_fix_calculation < number_of_prob_fix_calculations; each_prob_fix_calculation++){
        std::string Prob_fixation = choose_which_prob_fix_string(each_prob_fix_calculation);

// Initiate and Loop through Effective Population size
        for(int each_population_size = 0; each_population_size < number_of_pops; each_population_size++){
            long double N = pop_sizes[each_population_size];
            
// Initiate and Loop through Time               
            long double generation = 0.00;
            for(int each_generation_interval = 0; each_generation_interval < number_time_intervals; each_generation_interval++){
                long double time = generation;
                
//------Stochiometric Calculations----------------------------------
                long double concentration_hydrophobic_residues_per_tissue_in_balance = do_all_stoichiometric_calc_for_concentration_hp(k, a_total, b_total); 
                long double concentration_hydrophobic_residues_per_tissue_out_of_balance = do_all_stoichiometric_calc_for_concentration_hp(k, a_total/2, b_total); 

//------Fitness Calculations----------------------------------------
                std::vector <long double> fitness_values_based_on_hp_in_each_tissue;
                long double current_state = 0.0;
                for(int each_state = 0; each_state <= number_of_regulatory_regions; each_state++){
                    // Sum hp across all regulatory domains
                    long double sum_of_hp_per_state = sum_hp_across_tissues(number_of_regulatory_regions, current_state, concentration_hydrophobic_residues_per_tissue_in_balance, concentration_hydrophobic_residues_per_tissue_out_of_balance);
                    //Calculate the fitness of that state based on the [hp]
                    long double fitness_per_state = calculate_fitness_of_state(sum_of_hp_per_state);
                    fitness_values_based_on_hp_in_each_tissue.push_back(fitness_per_state);
//                    std::cout << "fitness: " << fitness_per_state << "\n\n";
                    current_state = current_state + 1;
                }

            //------Markov Matrix Calculations-----------------------------------                
                Eigen::MatrixXf probability_distribution_matrix(matrix_size, 1);
                probability_distribution_matrix = calculate_probability_distribution_across_states(fitness_values_based_on_hp_in_each_tissue, time, N, each_prob_fix_calculation);
                
                // Convert Eigen::Matrix to std::vector
                std::vector <long double> probability_distribution_vector(probability_distribution_matrix.data(), probability_distribution_matrix.data() + probability_distribution_matrix.size());
                        
            //------Print to .csv file--------------------------------------------                  
                myfile << Prob_fixation << ", "  << scalar << ",   "  << N << ",  "  << k << ",  "  <<time;
                for(int each_states_prob = 0; each_states_prob < matrix_size; each_states_prob++){
                    myfile << ",   " << probability_distribution_vector[each_states_prob];                    
                }
                myfile << "\n";

// Add to time interval                
                generation = generation + size_time_intervals;
            }
        }
    }

    myfile.close();
    return 0;
}

//------Function definitions------------------------------------------
std::string choose_which_prob_fix_string(int prob_fix_calculation){
//Tell loop whether to use prob of fix = 1/N or Sella and Hirsch equation    
    if(prob_fix_calculation == 0){
        std::string Prob_fixation_string = "Prop to [hp]";
        return Prob_fixation_string;
    }
    else if (prob_fix_calculation == 1){
        std::string Prob_fixation_string = "One over Ne";
        return Prob_fixation_string;
    }
    else{
        std::string Prob_fixation_string = "Houston, we have a problem";        
        return Prob_fixation_string;
    }
}

long double calculate_concentration_in_bound_form(long double keq, long double x_total, long double y_total){
//Calculate the concentration of of ab molecules in bound form using the quadradic equation
    long double a = keq ;
    long double b = -(keq*x_total+keq*y_total+1) ;
    long double c = keq*x_total*y_total;
    
    long double xy_1;
    xy_1 = (-b + sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);

    long double xy_2;
    xy_2 = (-b - sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
    
// choose which value to use for bound xy
    long double xy_bound; 
    if (xy_1 >= 0 && xy_2 > x_total && xy_2 > y_total){        
        xy_bound = xy_1;
    }
    else if (xy_2 >= 0 && xy_1 > x_total && xy_1 > y_total){
        xy_bound = xy_2;
    }
    else if (xy_2 >= 0 && xy_1 < 0){
        xy_bound = xy_2;
    }
    else if (xy_1 >= 0 && xy_2 < 0){
        xy_bound = xy_1;
    }
    else {
        std::cout << "ERROR";
    }
    return xy_bound;
}

long double calculate_concentration_x_free(long double x_total, long double x_bound){
//Calculate [X]free by subtracting [XY] from [X]tot        
    long double x_free_con = x_total - x_bound;
    return x_free_con;
}
    
long double calculate_concentration_hydrophobic_residues(long double x_free, long double y_free){
//Calculate [hp] by summing [X]free and [Y]free    
    long double hp =  x_free + y_free;
    return hp;
}

long double do_all_stoichiometric_calc_for_concentration_hp(long double keq, long double x_total, long double y_total){
//------Stochiometric Calculations----------------------------------
//Calculate the concentration of of xy molecules in bound form
    long double xy = calculate_concentration_in_bound_form(keq, x_total, y_total);
    
//Calculate [X]free and [Y]free by subtracting [XY] from [X]tot and [X]tot    
    long double x_free = calculate_concentration_x_free(x_total, xy);
    long double y_free = calculate_concentration_x_free(y_total, xy);
    
//Calculate [hp] by summing [A]free and [B]free    
    long double hp = calculate_concentration_hydrophobic_residues(x_free, y_free);
    return hp;
}

long double sum_hp_across_tissues(long double regulatory_regions, long double current_state, long double hp_per_tissue_in_balance, long double hp_per_tissue_out_of_balance){
//Sum hp across all regulatory domains
    long double unaffected_tissue = regulatory_regions-current_state;
    long double affected_tissue = current_state;
    long double sum_of_hp_across_tissues_for_state = unaffected_tissue*hp_per_tissue_in_balance + affected_tissue*hp_per_tissue_out_of_balance;
    return sum_of_hp_across_tissues_for_state;
}

long double calculate_fitness_of_state(long double sum_of_hp_across_tissues_for_state){
//Calculate the fitness of that state based on the [hp]
    long double fitness_of_state = 1/(scalar*(1+ sum_of_hp_across_tissues_for_state));
    return fitness_of_state;
}


long double which_prob_fix_calculation(long double fitness_current_state, long double fitness_next_state, long double N, int prob_of_fix_which){
// Decide which prob of fixation calculation to use (sella and hirsch or 1/N) and calculate it
    if(prob_of_fix_which == 0){
//------Probability of fixation proportional to hp---------------------  
        long double this_prob_fix = calculate_prob_fix_using_sella_hirsch(fitness_current_state, fitness_next_state, N);
        return this_prob_fix;
    }
    else if(prob_of_fix_which == 1){
//------Probability of fixation = 1/N-----------------------------------
        long double this_prob_fix = calculate_prob_fix_as_one_over_Ne(N); 
        return this_prob_fix;
    }
    else{
        std::cout << "ERROR";
        long double this_prob_fix = 0.00; 
        return this_prob_fix;        
    }
}


long double calculate_prob_fix_using_sella_hirsch(long double fitness_current_state, long double fitness_next_state, long double population_size){
//Calculate the Prob of fixation using the Sella and Hirsch equation 
    long double probability_of_fixation = (1-(fitness_current_state/fitness_next_state))/(1-(pow((fitness_current_state/fitness_next_state),population_size)))*population_size;
    return probability_of_fixation;
}

long double calculate_prob_fix_as_one_over_Ne(long double population_size){
//Calculate the Prob of fixation if you want to assume its 1/N
    long double probability_of_fixation = (1/population_size)*population_size;
    return probability_of_fixation;
}

Eigen::MatrixXf calculate_probability_distribution_across_states(std::vector <long double> fitness_of_states, long double time, long double N, int prob_fix_calculation){
//------Generator Rate Matrix (Q)--------------------------------------
//Initialize q matrix size
    Eigen::MatrixXf q_generator_rate_matrix(matrix_size, matrix_size);
    q_generator_rate_matrix << Eigen::MatrixXf::Zero(matrix_size, matrix_size);

// Edit empty q_generator_rate_matrix matrix with the proper formulas for the rate matrix
// Calculates the rate between each state by incorportating the probability of fixation from the Sella and Hirsch formula or 1/N
    int state_number;
    long double lost_regulatory_regions = 0.0;
    for (state_number = 0; state_number < number_of_regulatory_regions; state_number++){
        
        long double fitness_current_state = fitness_of_states[state_number]; 
        long double fitness_next_state = fitness_of_states[state_number+1];
        long double fjs = fitness_of_states[state_number+1];
        long double fjp = fitness_of_states[number_of_regulatory_regions];  
        
        long double prob_fix_transient = which_prob_fix_calculation(fitness_current_state, fitness_next_state, N, prob_fix_calculation);
        long double prob_fix_s = which_prob_fix_calculation(fitness_current_state, fjs, N, prob_fix_calculation);
        long double prob_fix_p = which_prob_fix_calculation(fitness_current_state, fjp, N, prob_fix_calculation);        
       
//----------------------------------------------------------------------
        
        if (state_number == 0){
            q_generator_rate_matrix(state_number, state_number) = -(prob_fix_transient*nucleotide_mutation_rate*nucleotide_length_regulatory_region*2*number_of_regulatory_regions_double + 2*prob_fix_p*nucleotide_mutation_rate*nucleotide_length_coding_region);
            q_generator_rate_matrix(state_number+1, state_number) = prob_fix_transient*nucleotide_mutation_rate*nucleotide_length_regulatory_region*2*number_of_regulatory_regions_double;
            q_generator_rate_matrix(number_of_regulatory_regions+1, state_number) = 2*prob_fix_p*nucleotide_mutation_rate*nucleotide_length_coding_region;

        }
        else if (state_number > 0 and state_number < number_of_regulatory_regions-1){
            q_generator_rate_matrix(state_number, state_number) = -(prob_fix_transient*nucleotide_mutation_rate*nucleotide_length_regulatory_region*(number_of_regulatory_regions_double-lost_regulatory_regions) + prob_fix_s*nucleotide_mutation_rate*nucleotide_length_regulatory_region*(number_of_regulatory_regions_double-lost_regulatory_regions) + prob_fix_p*nucleotide_mutation_rate*nucleotide_length_coding_region);
            q_generator_rate_matrix(state_number+1, state_number) = prob_fix_transient*nucleotide_mutation_rate*nucleotide_length_regulatory_region*(number_of_regulatory_regions_double-lost_regulatory_regions);
            q_generator_rate_matrix(number_of_regulatory_regions, state_number) = prob_fix_s*nucleotide_mutation_rate*nucleotide_length_regulatory_region*(number_of_regulatory_regions_double-lost_regulatory_regions);
            q_generator_rate_matrix(number_of_regulatory_regions+1, state_number) = prob_fix_p*nucleotide_mutation_rate*nucleotide_length_coding_region;

        }
        else if (state_number == number_of_regulatory_regions-1){
            q_generator_rate_matrix(state_number, state_number) = -(prob_fix_s*nucleotide_mutation_rate*nucleotide_length_regulatory_region + prob_fix_p*nucleotide_mutation_rate*nucleotide_length_coding_region + prob_fix_p*nucleotide_mutation_rate*nucleotide_length_regulatory_region);
            q_generator_rate_matrix(number_of_regulatory_regions, state_number) = prob_fix_s*nucleotide_mutation_rate*nucleotide_length_regulatory_region;
            q_generator_rate_matrix(number_of_regulatory_regions+1, state_number) = prob_fix_p*nucleotide_mutation_rate*nucleotide_length_coding_region + prob_fix_p*nucleotide_mutation_rate*nucleotide_length_regulatory_region; 
        }
        else{
            std::cout << "error";
        }
        lost_regulatory_regions = lost_regulatory_regions + 1;
    }  
// Prints q_generator_rate_matrix
//    std::cout << q_generator_rate_matrix << std::endl<< "\n";
    
//-------Probability (P) Matrix---------------------------------------
// Calculate q matrix * time and prints q*time matrix
    Eigen::MatrixXf qtime(matrix_size , matrix_size);
    qtime = q_generator_rate_matrix*time;
//    std::cout << qtime << std::endl<< "\n";

// Calculate the p matrix by exponentiation of the q*time matrix, P = e^(q*t)
    Eigen::MatrixXf probability_matrix(matrix_size , matrix_size);
    probability_matrix = qtime.exp();
//    std::cout << probability_matrix << "\n\n";
    
//-------calculate probability distribution across states--------------------------------------------------
// Initial state vector, with all 100 individuals starting as identical duplicates
    Eigen::MatrixXf initial_state_matrix(matrix_size, 1);
    initial_state_matrix << Eigen::MatrixXf::Zero(matrix_size,1);
    initial_state_matrix(0,0) = 100;
//    std::cout << initial_state_matrix << std::endl<< "\n";
//    std::vector <int> state = {100, 0, 0, 0, 0, 0};

// Calculate probability distribution from initial vector
    Eigen::MatrixXf probability_distribution(matrix_size, 1);
    probability_distribution = probability_matrix*initial_state_matrix;
    std::cout << probability_distribution << std::endl<< "\n";

    return probability_distribution.col(0);

}


//-------------------------------------------------------------------