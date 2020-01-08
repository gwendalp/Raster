#define BOOST_TEST_MODULE 

#include <boost/test/unit_test.hpp>
float area(double x1, double y1, double x2, double y2, double x3, double y3) 
    { 
        /*
          Paramètres entrée : les coordonnées des sommets d'un triangle
          
          Permet de retourner l'aire d'un triangle
        */
        return abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0); 
    } 


bool isInside(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y) 
    {   
        /* Calculate area of triangle ABC */
        float A = area (x1, y1, x2, y2, x3, y3); 
        
        /* Calculate area of triangle PBC */   
        float A1 = area (x, y, x2, y2, x3, y3); 
        
        /* Calculate area of triangle PAC */   
        float A2 = area (x1, y1, x, y, x3, y3); 
        
        /* Calculate area of triangle PAB */    
        float A3 = area (x1, y1, x2, y2, x, y); 
            
        /* Check if sum of A1, A2 and A3 is same as A */ 
        return (abs(A - (A1 + A2 + A3)) < 0.01); 
    } 



BOOST_AUTO_TEST_CASE(TestArea)
{
    
    double x1 = 0.; 
    double y1 = 0.;
    double x2 = 1.;
    double y2 = 0.;
    double x3 = 1.;
    double y3 = 1.;

    BOOST_CHECK_CLOSE(area(x1, y1, x2, y2, x3, y3), 0.5, 0.01);
}





