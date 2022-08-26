# Wilson_Liberles_2022
Dosage Balance as a time-dependent selective barrier to subfunctionalization

Files:
  -Program 1: sub_dos_2022_fig2.cpp
  -Program 2: sub_dos_2022_figures3-6.cpp
  -Program 3: sub_dos_2022_fig4b.cpp
  
  
Program 1: 
  -Calculates the concentration of exposed hydrophobic residues given values for keq, [A]total, [B]total
  -Creates data for figure 2
  -Inputs: (lines 32, 33; lines 64-66)
      -line 32: g < # kays values
      -line 33: h < # atots values
      -line 64: kays = {given keq values}
      -line 65: atots = {given [A]total}
      -line 66: btots = {given [B]total}

Program 2:
  -Calculates expected distribution of states for subfunctionalization-only model and subfunctionalization + Dosage model after whole-genome duplication
  -Currently written to produce data for figure 3 and 4a
	
	
  -Pseudo code
	
	
	//
        For each effective population size listed;
            For each generation until desired time point;
                Concentration of hydrophobic residues when in stoichiometric balance ([A]tot=[B]tot) = [A]free + [B]free (calculated from k, [A]total, [B]total);
                Concentration of hydrophobic residues when out of stoichiometric balance ([A]tot=2[B]tot) = [A]free + [B]free (calculated from k, [A]total, [B]total);
                For each state;
                    Unaffected tissues = z – current state;
                    Affected tissues  = current state;
                    Summation of hydrophobic residues across tissues = (unaffected tissues)*concentration hp when in stoichiometric balance + affected tissues * concentration hp when out of stoichiometric balance;
                    Fitness per state = 1/(scalar + sum of hp across tissues for state);
                Vector of fitness values for each state  = [fitness of state 0, fitness of state 1, etc];
                Probability Distribution Across States = using generator Markov matrix with vector of fitness values for each state, and initial state vector;
            
                
  ***To use this code need Eigen package***
  -Inputs: (lines 91 - 123)
  -AND change the initial state vector if need be (lines 376-380)
	
	
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
    
Program 3: 
  -Calculates expected distribution of states for subfunctionalization-only model and subfunctionalization + Dosage model after small-scale duplication. 
  -Currently written to produce data for figure 4b.
  -When calculating  sum of hp across tissues,first calculate [a]free and [b]free in balance and out of balance tissues and sum of [hp] = affected_tissue*hp_per_tissue_in_balance + unaffected_tissue*hp_per_tissue_out_of_balance;
  -and [A]tot and [B]tot is half that of WGD events...
  -long double a_total = 0.00000125;
  -long double b_total = 0.00000125; 

  ***To use this code need Eigen package***
  -Inputs: (lines 85 - 116)
  -AND change the initial state vector if need be (lines 368-373)
	
	
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
