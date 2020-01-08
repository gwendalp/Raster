#include "delaunator.hpp"
#include <cstdio>
#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <map>
#include "WGS84toCartesian.hpp"
#include "colourmanager.h"
#include <algorithm>
using namespace std;

/* Variable globale utile */
vector<double> coords;
vector<double> X;
vector<double> Y;
vector<double> Z;
map<tuple<double, double>, double> xy2z;
map<int, vector<int>> label;


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
int nb_pixel_x;
int nb_pixel_y;



void printProgress (double percentage) //Affichage du barre de progrès dans le terminal
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

void Mercator(double &lx, double &ly) //Projection Mercator
    { 
        /*
          Entrées/Sorties : lx qui est la coordonnées X du fichier
                            ly qui est la coordonnées Y du fichier
        */
        float lx_0 = 48.29762887;
        float ly_0= -004.41737340;

        array<double,2> pos{wgs84::toCartesian({lx_0,ly_0},{lx,ly})};
        ly = pos[0]; // Récupère la coordonnée en x
        lx = pos[1]; // Récupère la coordonnée en y
    }

void File_read(char *file_name,vector<double> &X, vector<double> &Y, vector<double> &Z)
    {   
        /*
          Paramètres entrée : file_name qui permet de récuperer du terminal le nom du fichier à ouvrir
                              X,Y,Z les coordonnées du fichier placées en référence
        */
       //Déclaration des variables locales
        double minx;
        double miny;
        double lx, ly, lz;

        ifstream Flux(file_name);
        int c = 0;
        if(Flux.good())
        {
            while(Flux.good())
            {
                Flux >> lx;
                Flux >> ly;
                Flux >> lz;

                if (c % 10 == 0)
                {
                    Mercator(lx,ly);
                    X.push_back(lx);
                    Y.push_back(ly);
                    Z.push_back(lz);
                }
                c+=1;
            }       
        }

        minx = *min_element(X.begin(), X.end());
        miny = *min_element(Y.begin(), Y.end());
        
        for (int i = 0; i < X.size(); i++) 
        {
            //permet de placer un offset pour mettre le minimum de X et Y à 0
            X[i] -= minx;
            Y[i] -= miny;
        }        
    }


void Labelize(vector<vector<double>> &triangle, map<int, vector<int>> &label, double pas_zonex, double pas_zoney)
    {
        /*
          Paramètres entrée : Triangle est un vector contenant l'élément triangle qui est un vector
                              Map est le dictionnaire, la clé est le numéro de maille du maillage
                                                       la valeur est un vector contenant l'indice du triangles dans le vector Triangle
        */
        for(int i=0; i<triangle.size(); i++)
        {
            int zone;
            /*
              Zone permet de désigner la zone du maillage correspondant à un des sommets du triangle
              La condition if permet d'éviter les doubles d'indices de triangles pour le maillage
            */
            zone = int(int(triangle[i][0]/pas_zonex)*10 + int(triangle[i][1]/pas_zoney));
            if(find(label[zone].begin(), label[zone].end(), i) != label[zone].end()) {} 
            else 
            {
                label[zone].push_back(i);
            }
            
            
            zone = int(int(triangle[i][2]/pas_zonex)*10 + int(triangle[i][3]/pas_zoney));
            if(find(label[zone].begin(), label[zone].end(), i) != label[zone].end()) {} 
            else 
            {
                label[zone].push_back(i);
            }


            zone = int(int(triangle[i][4]/pas_zonex)*10 + int(triangle[i][5]/pas_zoney));
            if(find(label[zone].begin(), label[zone].end(), i) != label[zone].end()) {} 
            else 
            {
                label[zone].push_back(i);
            }
        }
    }
void post_traitement(int npx, int npy, vector< vector<double>> &img)
    {
        /*
          Cette fonction permet de traiter après coup les petites erreurs dûes au maillage
        */
        int compteur = 0;
        for (int i = 1; i < npx - 2;  i++)
        {
            for (int j = 1; j < npy -2; j++)
            {
                if(img[i][j] == 0)
                {   
                    compteur += 1;
                    int c = 0;
                    if(img[i-1][j] == 0) {c += 1;}
                    if(img[i-1][j+1] == 0) {c += 1;}
                    if(img[i-1][j-1] == 0) {c += 1;}
                    if(img[i][j+1] == 0) {c += 1;}
                    if(img[i][j-1] == 0) {c += 1;}
                    if(img[i+1][j+1] == 0) {c += 1;}
                    if(img[i+1][j-1] == 0) {c += 1;}
                    if(img[i+1][j] == 0) {c += 1;}

                    

                    if(c < 5) // Sur un ensemble convexe : augmenter cette valeur, la diminuer pour un concave.
                    {
                        img[i][j] = max(img[i][j+1],max(img[i+1][j],max(img[i-1][j+1],max(img[i+1][j-1],max(img[i+1][j+1],max(img[i-1][j],max(img[i-1][j-1],img[i][j-1])))))));
                    }

                    else 
                    {
                        
                        img[i][j] = 0;
                    }
                }        
            }
        }

        //cout << compteur << endl;    
    }

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
void giveDepth(double px, double py,vector<vector<double>> &triangle, map<int, vector<int>> &label, double pas_zonex, double pas_zoney, double hx, vector<vector <double>> &img)
    {
        /*
            Paramètres entrée : px et py les coordonnées du pixel
                                triangle qui contient tous les coordonnées de chaque triangle
                                map étant notre dictionnaire reliant le maillage à l'indice des triangles dans "triangle"
                                img est la matrice des profondeurs
        */
        double x = px*hx;
        double y = py*hx;
        int npx = img[0].size();
        int npy = img[1].size();
        int zone = int(int(x/pas_zonex)*10 + int(y/pas_zoney));

        for(int i = 0; i<label[zone].size(); i++)
        {
            if (isInside(triangle[label[zone][i]][0], triangle[label[zone][i]][1], triangle[label[zone][i]][2], triangle[label[zone][i]][3], triangle[label[zone][i]][4], triangle[label[zone][i]][5], x, y))
            {
                double z1 = xy2z[{triangle[label[zone][i]][0], triangle[label[zone][i]][1]}];
                double z2 = xy2z[{triangle[label[zone][i]][2], triangle[label[zone][i]][3]}];
                double z3 = xy2z[{triangle[label[zone][i]][4], triangle[label[zone][i]][5]}];  
                img[px][py] = 0.33333*(z1 + z2 + z3) ;

                return; 
            }
        }
    }


int main(int argc, char * argv[]) 
{
    char *file_name;
    file_name = argv[1];
    File_read(file_name,X, Y, Z);

    for (int i = 0; i < X.size(); i++)
    {      
        /*
          Crée le dictionnaire qui associe à chaque couple (X,Y) une profondeur Z
        */ 
        coords.push_back(X[i]);
        coords.push_back(Y[i]);
        xy2z[{X[i], Y[i]}] = Z[i];
    }

    double minz = *min_element(Z.begin(), Z.end());
    double minx = *min_element(X.begin(), X.end());
    double miny = *min_element(Y.begin(), Y.end());
    double maxz = *max_element(Z.begin(), Z.end());
    double maxx = *max_element(X.begin(), X.end());
    double maxy = *max_element(Y.begin(), Y.end());
  
    /*
      Bibliothèque permettant de créer une color map : à chaque profondeur on associe un tuple RGB
    */
    ColourManager::Init_ColourManager(); 

    ColourManager manager(minz, maxz);
    ColourMap cmap("cmap");             // Couleur de la légende :
    cmap.addColour(142, 68, 173,1.0);   // Wisteria
    cmap.addColour(41, 128, 185, 1.0);  // Belize Hole
    cmap.addColour(46, 204, 113, 1.0);  // Nephritis
    cmap.addColour(241, 196, 15, 1.0);  // Sun Flower
    cmap.addColour(230, 126, 34, 1.0);  // Carrot
    cmap.addColour(235, 59, 90,1.0);    // Crimson 

    // Setting the ColourMap
    ColourManager::setCurrentColourMap(cmap);
    Colour c = manager.getInterpolatedColour(0);
    c = manager.getInterpolatedColour(0);

    /*
    Récupère le nombre de pixel souhaité
    */
    nb_pixel_x = atoi(argv[2]); //Permet de convertir un char* en int
    double hx = (maxx)/nb_pixel_x;
    double hy = hx; 
    int nb_pixel_y = int(maxy/hy);


    vector<vector <double>> img;

    for (int i = 0; i < nb_pixel_x ; i++)
    {
        /*
          img est la matrice des profondeurs
          Ces boucles permettent d'initialiser cette matrice img
        */
        img.push_back({});
        for (int j = 0; j < nb_pixel_y; j++)
        {
            img[i].push_back(0);
        }   
    }
    
    /*
      Triangulation of Delaunator happens here
      Récupère les coordonnées des sommets du triangle afin de l'ajouter en tant que vector dans le vector de vector triangles
    */
    delaunator::Delaunator d(coords);

    vector<vector<double>> triangle;
     for(std::size_t i = 0; i < d.triangles.size(); i+=3) 
    {
            double x1 = d.coords[2 * d.triangles[i]];       
            double y1 = d.coords[2 * d.triangles[i] + 1];   
            double x2 = d.coords[2 * d.triangles[i + 1]];   
            double y2 = d.coords[2 * d.triangles[i + 1] + 1];
            double x3 = d.coords[2 * d.triangles[i + 2]];    
            double y3 =  d.coords[2 * d.triangles[i + 2] + 1]; 
            
            
            
            if((std::max(abs(x1-x2), abs(x2 -x3)) < 0.1*maxx) or (std::max(abs(y1-y2), abs(y2 -y3)) < 0.1*maxy))
            {
            triangle.push_back({x1, y1, x2, y2, x3, y3});
            }
    }
    
    //pas pour le maillage
    double paszonex = int(maxx/100);
    double paszoney = int(maxy/100);

    //Appel de la fonction permettant le maillage
    Labelize(triangle, label, paszonex, paszoney);
        
    for (int i = 1; i < nb_pixel_x; i++)
    {
        for (int j = 1; j < nb_pixel_y; j++)
        {
            /*
            Permet de retourner la profondeur associée au pixel
            */
            giveDepth(i, j, triangle, label, paszonex, paszoney,  hx, img);

            


        }

        printProgress(i*float(1./float(nb_pixel_x)));
    }
    
    

    // Permet de lisser, les pixels non attribués 

    post_traitement(nb_pixel_x, nb_pixel_y, img);
    post_traitement(nb_pixel_x, nb_pixel_y, img);
    post_traitement(nb_pixel_x, nb_pixel_y, img);
    
    
    
    
    /*
      Permet de remplir le fichier .ppm pour afficher l'image
    */
    ofstream image;
    image.open("image.ppm", ofstream::out | ofstream::trunc);
    image << "P6" << endl << nb_pixel_y << " " << nb_pixel_x << endl << 255 << endl;
    for (int i = nb_pixel_x-1; i >-1; i--)
    {
        for (int j = 0; j <nb_pixel_y;j++)
        {
            if (img[i][j] == 0)
            {
                image << char(0) << char(0) <<  char(0);
            }
            else
            {
                c = manager.getInterpolatedColour(img[i][j]);
                image <<  char(c.getIntR()) << char(c.getIntG()) << char(c.getIntB());
            }
        }
        
    }
}