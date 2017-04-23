#include "XSec.h"
#include "Solver.h" // Linterpolate
#include <iostream>

double lookup_XSec:: xs( const double E ) const
{ 
    double E1,E2,XS1,XS2,XS;
    for (int i = 0; i<Edata.size(); i++)
    {
        if (Edata[i]>=E )
        {
            E1=Edata[i]; E2=Edata[i+1];
            XS1=XSdata[i];XS2=XSdata[i+1];
            break;
        }      
    }
    XS=Linterpolate(E,E1,E2,XS1,XS2);
    //std::cout<<"E= "<<E<<" XS= " <<XS<<std::endl;
    //std::cout<<"E1= "<<E1<<" E2= " <<E2<<std::endl;
    return XS; 
}