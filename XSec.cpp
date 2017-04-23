#include "XSec.h"
#include "Solver.h" // Linterpolate
#include <iostream>

double lookup_XSec:: xs( const double E ) const
{ 
    double E1,E2,XS1,XS2,XS;
    int bin = Binary_Search(E,Edata);
    
    E1=Edata[bin]; E2=Edata[bin+1];
    XS1=XSdata[bin];XS2=XSdata[bin+1];
            
    XS=Linterpolate(E,E1,E2,XS1,XS2);
    //std::cout<<"E= "<<E<<" XS= " <<XS<<std::endl;
    //std::cout<<"E1= "<<E1<<" E2= " <<E2<<std::endl;
    return XS; 
}