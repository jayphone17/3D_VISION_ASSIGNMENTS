#ifndef CLOUD_SEGMENTATION_H
#define CLOUD_SEGMENTATION_H

#include "types.h"
#include "Ply.hpp"
#include "Visualization_Tools.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

using I::misc::Ply;
typedef unsigned char byte;
const double PI = 3.141592;





class Cloud_segmentation
{
	
	struct HPoint {
		Point_d position;
		Vector_3 normal;
		int primitive_index;
	};


 

public:

	vector<HPoint> HPS;
	std::vector < std::vector < int > > spherical_neighborhood;
	std::vector < std::vector < int > > plane_point_index;
	std::vector < Plane_3 > extracted_planes;
	bounding_3 BBox_scan;
	double BBox_diagonal;
	



bool load_points(string filename)
	{
		cout<<"Loading ";
		cout <<filename;

		HPS.clear();
		Ply ply;

		if (!ply.open(filename)){ return false;}

		for (Ply::ElementsIterator it = ply.elements_begin(); it != ply.elements_end(); ++it){
			const Ply::Element& element = *it;

			if (element.name() != "vertex"){
				if (!ply.skip(element)){ ply.close(); return false;}
				continue;
			}

			/* Check the properties are what we expect */
			if ((element.property(0).name() != "x")||(element.property(1).name() != "y")||(element.property(2).name() != "z")){cerr << "Unexpected vertex properties in the PLY file" << endl; ply.close(); return false;}


			size_t num_vertices = element.count();
			HPS.resize(num_vertices);

			

			if (element.num_properties() == 6) {
				for (size_t i=0; i<num_vertices; i++){
					double x,y,z;
					double nx,ny,nz;

					if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){cerr << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
					if ((!ply.read(element.property(3), nx))||(!ply.read(element.property(4), ny))||(!ply.read(element.property(5), nz))){cerr << "error while reading attribut " << i+1 << endl; ply.close(); return false; }

					HPS[i].position = Point_d(x,y,z);
					HPS[i].normal = Vector_3(nx,ny,nz);
				}
			}

			

			else {cerr << "Unexpected vertex properties in the PLY file" << endl; ply.close(); return false;}
		}

		ply.close();

		cout<<" ("<<HPS.size()<<" points)"<<endl;   
		return true;
	}






bool Compute_Knearest_neighbors(int K)
	{

		cout<<"Computing the K-nearest neighbors";

		//1-Neighborhood computation and reset the attributes of the structure points
		std::map<Point_d,int> map_indice_point;

		for(int i=0;i<HPS.size();i++){
			std::vector < int > plane_index_list_tmp;
			Point_d pt=HPS[i].position;
			map_indice_point[pt]=i;
		}

		std::list<Point_d> list_points;
		for(int i=0;i<(int)HPS.size();i++){
			Point_d pt=HPS[i].position;
			list_points.push_back(pt);
		}

		Tree tree(list_points.begin(), list_points.end());

		for(int i=0;i<(int)HPS.size();i++){
			Point_d query=HPS[i].position;
			Neighbor_search search(tree, query, K+1);

			std::vector < int > index_of_neighbors;
			for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
			//	if(std::sqrt(it->second)<RADIUS_NEIGHBORS_MAX){
					std::map<Point_d,int>::iterator iter =map_indice_point.begin();
					iter= map_indice_point.find(it->first);
					if( iter != map_indice_point.end() && iter->second!=i ) index_of_neighbors.push_back(iter->second);
			//	}else{break;}
			}
			spherical_neighborhood.push_back(index_of_neighbors);
		}

		//2-creation of the bounding box
		BBox_scan = CGAL::bounding_box(list_points.begin(), list_points.end());
		BBox_diagonal=sqrt( pow(BBox_scan.xmax()-BBox_scan.xmin(),2) +  pow(BBox_scan.ymax()-BBox_scan.ymin(),2) +  pow(BBox_scan.zmax()-BBox_scan.zmin(),2) );
	
		cout<<endl<<endl;

		return true;
	}




	


bool region_growing(float epsilon, int Nmin){
			

		cout<<"Region growing ";


		//Initialization structures
		plane_point_index.clear();
		for(int i=0; i<(int)HPS.size();i++) HPS[i].primitive_index=-1;
		int class_index=-1;
	

		// for each point, if not inliers yet..
		for(int i=0; i<(int)HPS.size();i++){

			if(HPS[i].primitive_index==-1){

				//update the index of primitive
				class_index++;
				HPS[i].primitive_index=class_index;

				//characteristics of the seed
				Vector_3 normal_seed=HPS[i].normal;
				Point_d pt_seed=HPS[i].position;
				Plane_3 optimal_plane(pt_seed,normal_seed);

				//initialization containers
				std::vector < int > index_container; index_container.push_back(i);
				std::vector < int > index_container_former_ring; index_container_former_ring.push_back(i);
				std::list < int > index_container_current_ring;



				//propagation
				bool propagation=true;
				do{

					propagation=false;

					// TO START FROM HERE (YOU)
					//-----------------------
					
 

					//-----------------------

					// update containers
					index_container_former_ring.clear();
					for(std::list < int >::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it){
						index_container_former_ring.push_back(*it);
						index_container.push_back(*it);
					}
					index_container_current_ring.clear();

				}while(propagation);
				

				//Test the number of inliers -> reject if inferior to Nmin
				if(index_container.size()>=Nmin) plane_point_index.push_back(index_container);
				else{ 
					class_index--;
					HPS[i].primitive_index=-1;
					for(int k=0;k<(int)index_container.size();k++) HPS[index_container[k]].primitive_index=-1; 
				}
			} 
		}
		
		cout<<endl<<"-> "<<plane_point_index.size()<<" primitives"<<endl<<endl;
		return true;
	}










protected: 


};


#endif 
