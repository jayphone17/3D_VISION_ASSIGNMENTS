
#include <iostream>
#include <vector>
#include <cstring>


#include "types.h"
#include "cloud_segmentation.h"





int main(int argc, char **argv)
{
	

// -------  Loading data and structure initialization -------//
// ----------------------------------------------------------//
Cloud_segmentation C;

//points and normals
if(!C.load_points(argv[1])) return 0;


/*
//maximal distance error
float epsilon2=atof(argv[2]);  

//minimal number of inliers
int Nmin2=atoi(argv[3]);

*/

//Compute K-nearest neighbors 
C.Compute_Knearest_neighbors(20);


// ----------------------------------------------------------//
// ----------------------------------------------------------//







// --------------------- Region Growing ---------------------//
// ----------------------------------------------------------//
bool redo=true;

while(redo){

	double epsilon;
	double Nmin;
	cout<<endl<<"Give value for epsilon: ";
	std::cin>>epsilon;
	cout<<"Give minimal number of inliers: ";
	std::cin>>Nmin;
	cout<<endl;


	//Region Growing HERE
	C.region_growing(epsilon,Nmin);

	//saving geometric primitives
	if(C.plane_point_index.size()>0) save_envelops(C);

	cout<<endl<<"Relaunch ? (yes=1, no=0): ";
	std::cin>>redo;
}


// ----------------------------------------------------------//
// ----------------------------------------------------------//




cout << endl<< "END" << endl;

	return 0;
}