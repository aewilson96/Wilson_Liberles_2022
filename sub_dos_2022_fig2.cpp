/*
Created by Amanda Erin Wilson on 03/18/2022
Published to Github 08/26/2022

g++ -I eigen/ sub_dos_2022_fig2.cpp -o program1

Calculates the concentration of exposed hydrophobic residues given values for keq, [A]total, [B]total
Creates data for figure 2

*/

#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include <fstream>

// function declaration
long double quad(long double k, long double atot, long double btot);

std::tuple <long double, long double, long double> concentrations(long double k, long double atot, long double btot, long double ab);

std::tuple <long double, long double, long double> katabat_list(int ik, int ia, int ib);

long double ratio_one(long double a, long double b);

int main() {
    std::ofstream myfile;
    myfile.open ("dimer_hp_sim_output_twodimensional2.csv");
    myfile << "keq,  Atot,  Btot,    Afree,   Bfree,    hp,    AB,  Imbalance,  log(keq),   log(hp)\n";
    
    for (int g = 0; g < 5; g++){
        for(int h = 0; h < 55; h++){
                std::tuple <long double, long double, long double> katabat = katabat_list(g, h, 0);

            // initialize number keq [Atot] [Btot]
                long double k = std::get<0>(katabat); // keq
                long double at = std::get<1>(katabat); // [Atot]
                long double bt = std::get<2>(katabat); // [Btot]                   
                
                long double z = quad(k, at, bt);
                
//                    long double ratio = at/bt;
                long double ratio = (1- ratio_one(at, bt))*100;
                
                
                std::tuple <long double, long double, long double> cs = concentrations(k, at, bt, z);
                
                long double logkeq = log10(k);
                long double hp = std::get<2>(cs);
                long double loghp = log10(hp);
                                    
                myfile << k << ",    "  << at << ",    "  << bt << ",    "  << std::get<0>(cs) << ",    "  << std::get<1>(cs) << ",    "  << std::get<2>(cs) << ",    "  <<z<< ",    " << ratio << ",    "<< logkeq << ",    "<< loghp <<"\n";
        }     
    }
    myfile.close();
            return 0;
}

// function definitions
//Give you k, [Atot], and [Btot] from lists
std::tuple <long double, long double, long double> katabat_list(int ik, int ia, int ib){
//  4, 55, 1 (output_twodimensional2)
    std::vector <long double> kays = {1000000.0, 100000000.0, 10000000000.0,  1000000000000.0, 100000000000000.0};
    std::vector <long double> atots = {0.0001, 0.0002, 0.0003,  0.0004, 0.0005, 0.0006, 0.0007, 0.0008,  0.0009, 0.00001, 0.00002,  0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009, 0.000001, 0.000002, 0.000003, 0.000004, 0.000005,  0.000006, 0.000007, 0.000008, 0.000009, 0.0000001, 0.0000002, 0.0000003, 0.0000004, 0.0000005, 0.0000006, 0.0000007, 0.0000008, 0.0000009, 0.00000001, 0.00000002, 0.00000003, 0.00000004, 0.00000005, 0.00000006, 0.00000007, 0.00000008, 0.00000009, 0.000000001, 0.000000002, 0.000000003, 0.000000004, 0.000000005,  0.000000006, 0.000000007,  0.000000008,  0.000000009, 0.0000000001};
    std::vector <long double> btots = {0.0000001};    
            
    long double k1 = kays[ik];
    long double at1 = atots[ia];
    long double bt1 = btots[ib];
    
    std::tuple <long double, long double, long double> katabat1 = {k1, at1, bt1};
    return katabat1;
}

// quadradic equations to calculate z/zed/possible values of [AB]
long double quad(long double k, long double atot, long double btot){
    //initialize 
    long double a = k ;
    long double b = -(k*atot+k*btot+1) ;
    long double c = k*atot*btot;
                    
    long double zee;
    zee = (-b + sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
//    std::cout << zee << std::endl;

    long double zed;
    zed = (-b - sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
//    std::cout << zed << std::endl;
    
// choose which z to use
    long double zs; 
    if (zee >= 0 && zed > atot && zed > btot){        
        zs = zee;
    }
    else if (zed >= 0 && zee > atot && zee > btot){
        zs = zed;
    }
    else if (zed >= 0 && zee < 0){
        zs = zed;
    }
    else if (zee >= 0 && zed < 0){
        zs = zee;
    }
    else {
        std::cout << "ERROR";
    }
    return zs;
    }


// Calculate x [A] and y [B]
std::tuple <long double, long double, long double> concentrations(long double k, long double atot, long double btot, long double ab){
// calculate [A]free and [B]free by subtracting [AB] from [A]tot and [B]tot    
    long double a_free_con = atot - ab;
    
    long double b_free_con = btot - ab;
    
// calculate [hp] by summing [A]free and [B]free    
    long double hp1 =  a_free_con + b_free_con;
//    std::cout << hp1 << std::endl;
    
    std::tuple <long double, long double, long double> cons = {a_free_con, b_free_con, hp1};
    return cons;
}

long double ratio_one(long double a, long double b){
    long double ratio;
    if (a >= b){
        ratio = b/a;
    }
    else if(b > a){
        ratio = a/b;
    }
    else{
        ratio = 0;
    }
    return ratio;
}