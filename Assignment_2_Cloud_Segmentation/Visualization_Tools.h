#ifndef VISUALIZATION_TOOLS_H
#define VISUALIZATION_TOOLS_H

#include "types.h"
#include "cloud_segmentation.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include <string>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/IO/print_wavefront.h>

#include <boost/filesystem.hpp>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <cmath>
#include <math.h>
#include <map>

#define RAD2DEG		180.0/M_PI
#define DEG2RAD		M_PI/180.0
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef Kernel::Point_2 Point_2;



#pragma warning(default:4244 4103 4267 4396)

using std::vector;

#define PARAM_OFSTREAM_TEXTE std::ios::out
#define OFSTREAM_TEXTE(nomvar, nomfic) std::ofstream nomvar (nomfic, PARAM_OFSTREAM_TEXTE) ; nomvar.setf(std::ios::fixed);



bool colorplanset2ply( char * filename, std::vector < std::vector < Point_d > > & pt, std::vector< CGAL::Color > & colors)
	{

	std::cout <<"Saving planes to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
	
	int nb_pts=0;
	for(int i=0;i<pt.size();i++){nb_pts=nb_pts+pt[i].size();}
		

	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<nb_pts<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property uchar red"<<std::endl;
		fic<<"property uchar green"<<std::endl;
		fic<<"property uchar blue"<<std::endl;
		fic<<"element face "<<pt.size()<<std::endl;
		fic<<"property list uchar int vertex_index"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<pt.size();i++){
			for(int j=0;j<pt[i].size();j++){
				fic<<pt[i][j].x()<<" "<<pt[i][j].y()<<" "<<pt[i][j].z()<<" "<<(int)colors[i].r()<<" "<<(int)colors[i].g()<<" "<<(int)colors[i].b()<<std::endl;
			}
		}
		int cont=0;
		for(int i=0;i<pt.size();i++){
			fic<<pt[i].size()<<" ";
			for(int k=cont;k<cont+pt[i].size();k++){
				fic<<k<<" ";}
		fic<<std::endl;
		cont=cont+pt[i].size();
		}

fic<<std::endl;	

return true;
	}


bool color_inliers( char * filename, std::vector < std::vector < Point_d > > & pt, std::vector< CGAL::Color > & colors)
	{

	std::cout <<"Saving planes to "<<filename<<std::endl;
	OFSTREAM_TEXTE(fic,filename);
	
	int nb_pts=0;
	for(int i=0;i<pt.size();i++){nb_pts=nb_pts+pt[i].size();}
		

	//header
		fic<<"ply"<<std::endl;
		fic<<"format ascii 1.0"<<std::endl;
		fic<<"comment author: F. Lafarge"<<std::endl;
		fic<<"element vertex "<<nb_pts<<std::endl;
		fic<<"property float x"<<std::endl;
		fic<<"property float y"<<std::endl;
		fic<<"property float z"<<std::endl;
		fic<<"property uchar red"<<std::endl;
		fic<<"property uchar green"<<std::endl;
		fic<<"property uchar blue"<<std::endl;
		fic<<"end_header"<<std::endl;

		//vertex list

		for(int i=0;i<pt.size();i++){
			for(int j=0;j<pt[i].size();j++){
				fic<<pt[i][j].x()<<" "<<pt[i][j].y()<<" "<<pt[i][j].z()<<" "<<(int)colors[i].r()<<" "<<(int)colors[i].g()<<" "<<(int)colors[i].b()<<std::endl;
			}
		}

fic<<std::endl;	

return true;
	}


template <class Cloud_segmentation>
bool save_envelops(Cloud_segmentation C)
	{
		std::vector < std::vector < Point_d > > point3dres_graham;
		std::vector < std::vector < Point_d > > point3dres_alpha_shape;
		
		std::vector < CGAL::Color > color_col_graham;
		std::vector < CGAL::Color > color_col_alpha_shape;

		int val_lab=0;
		int counting1=0;
		C.extracted_planes.clear();

		std::vector<CGAL::Color> vector_of_colors;
		for(int i=0; i<(int)C.plane_point_index.size();i++){
			int red=(int)floor((double)126*rand()/RAND_MAX)+130;
			int green=(int)floor((double)126*rand()/ RAND_MAX)+130;
			int blue=(int)floor((double)126*rand()/ RAND_MAX)+130;
			CGAL::Color col(std::max(0,std::min(255,red)),std::max(0,std::min(255,green)),std::max(0,std::min(255,blue)),120);
			vector_of_colors.push_back(col);
		}
		
		std::vector < std::vector < Point_d > > listpts;

		for(int i=0; i<(int)C.plane_point_index.size();i++){
			
			std::vector < Point_d > listp;
			for(int k=0; k<(int)C.plane_point_index[i].size(); k++){
				int yy=C.plane_point_index[i][k];
				Point_d pt=C.HPS[yy].position;
				listp.push_back(pt);
			}
			listpts.push_back(listp);

			val_lab++;
			Plane_3 plane;
			double erro=linear_least_squares_fitting_3(listp.begin(),listp.end(),plane, CGAL::Dimension_tag<0>());
			C.extracted_planes.push_back(plane);

			CGAL::Color col=vector_of_colors[i];
			
			
			Vector_3 vortho=plane.orthogonal_vector();
			Vector_3 b1=plane.base1();
			Vector_3 b2=plane.base2();
			FT norme1=sqrt(pow(b1.x(),2)+pow(b1.y(),2)+pow(b1.z(),2)); if(norme1<1e-7){norme1=1e-7;}
			FT norme2=sqrt(pow(b2.x(),2)+pow(b2.y(),2)+pow(b2.z(),2)); if(norme2<1e-7){norme2=1e-7;}

			Vector_3 vb1=(1/norme1)*b1;
			Vector_3 vb2=(1/norme2)*b2;

			std::vector< Point_d > poinn;
			std::vector< Point_2d > poin2;

			for(int k=0; k<(int)C.plane_point_index[i].size();k++){
				int ind=C.plane_point_index[i][k];
				Point_d pt=C.HPS[ind].position;
				Point_d ptp=plane.projection(pt);
				FT X1=vb1.x()*ptp.x()+vb1.y()*ptp.y()+vb1.z()*ptp.z();
				FT X2=vb2.x()*ptp.x()+vb2.y()*ptp.y()+vb2.z()*ptp.z();
				Point_2d ptp2(X1,X2);
				poinn.push_back(ptp);
				poin2.push_back(ptp2);
			}

			std::vector < Point_2d > poin2res_graham;
			std::vector < Point_2d > poin2res_alphashape;
			CGAL::ch_graham_andrew(poin2.begin(),poin2.end(),std::back_inserter(poin2res_graham));
			Alpha_Shape AlphaShape(poin2.begin(),poin2.end(),FT(0.005*C.BBox_diagonal),Alpha_Shape::REGULARIZED);   
			
			std::vector < std::pair < Point_2d, Point_2d> > pairs_edges;
			std::vector < bool > active_pairs;
			for (Alpha_Shape::Alpha_shape_edges_iterator it = AlphaShape.alpha_shape_edges_begin(); it != AlphaShape.alpha_shape_edges_end(); ++it){
				Alpha_Shape::Edge edg= *it;
				Point_2d pe1=edg.first->vertex(edg.first->ccw(edg.second))->point();
				Point_2d pe2=edg.first->vertex(edg.first->cw(edg.second))->point();
			
				std::pair < Point_2d, Point_2d> pair_ed;
				pair_ed.first=pe1;
				pair_ed.second=pe2;
				pairs_edges.push_back(pair_ed);
				active_pairs.push_back(true);
			}

			std::vector < Point_2d > ranked_pts;
			if(pairs_edges.size()>0){
				ranked_pts.push_back(pairs_edges[0].first); 
				ranked_pts.push_back(pairs_edges[0].second); 
				active_pairs[0]=false;
				Point_2d pt_to_find=pairs_edges[0].second;
				do{
					for(int iter=1;iter<pairs_edges.size();iter++){
						if(active_pairs[iter]){
							if(pt_to_find == pairs_edges[iter].first) { ranked_pts.push_back( pairs_edges[iter].second);  pt_to_find=pairs_edges[iter].second; active_pairs[iter]=false; }
							else if(pt_to_find == pairs_edges[iter].second) { ranked_pts.push_back( pairs_edges[iter].first); pt_to_find=pairs_edges[iter].first; active_pairs[iter]=false; }
						}
					}
				}while(pt_to_find!=pairs_edges[0].first);
			}

			Polygon_2 poly2(ranked_pts.begin(),ranked_pts.end());
			CDT2 cdt;

			for(Polygon_2::Vertex_iterator it2 = poly2.vertices_begin () ; it2 != poly2.vertices_end ();++it2) cdt.insert(*it2);

			CGAL::make_conforming_Delaunay_2(cdt);

			std::vector < Polygon_2 >  partition_polys;
			for ( CDT2::Finite_faces_iterator  vit = cdt.finite_faces_begin(); vit != cdt.finite_faces_end(); ++vit){        
                    vector<Point_2d> p(3);
                    p[0] = vit->vertex (0)->point();
                    p[1] = vit->vertex (1)->point();
                    p[2] = vit->vertex (2)->point();
                    Point_2d bary( (p[0].x()+ p[1].x() + p[2].x())/3.,(p[0].y()+ p[1].y() + p[2].y())/3.);
    
					if (CGAL::bounded_side_2(poly2.vertices_begin (), poly2.vertices_end (),bary, Kernel()) != CGAL::ON_UNBOUNDED_SIDE)
                    {
                        Polygon_2 po(p.begin(),p.end());
                        partition_polys.push_back(po);
						vector<Point_d> pd(3);
						FT X1=p[0].x()*vb1.x()+p[0].y()*vb2.x()-plane.d()*vortho.x();
						FT Y1=p[0].x()*vb1.y()+p[0].y()*vb2.y()-plane.d()*vortho.y();
						FT Z1=p[0].x()*vb1.z()+p[0].y()*vb2.z()-plane.d()*vortho.z();
						pd[0] = Point_d(X1,Y1,Z1); 
						FT X2=p[1].x()*vb1.x()+p[1].y()*vb2.x()-plane.d()*vortho.x();
						FT Y2=p[1].x()*vb1.y()+p[1].y()*vb2.y()-plane.d()*vortho.y();
						FT Z2=p[1].x()*vb1.z()+p[1].y()*vb2.z()-plane.d()*vortho.z();
						pd[1] = Point_d(X2,Y2,Z2); 
						FT X3=p[2].x()*vb1.x()+p[2].y()*vb2.x()-plane.d()*vortho.x();
						FT Y3=p[2].x()*vb1.y()+p[2].y()*vb2.y()-plane.d()*vortho.y();
						FT Z3=p[2].x()*vb1.z()+p[2].y()*vb2.z()-plane.d()*vortho.z();
						pd[2] = Point_d(X3,Y3,Z3); 
						point3dres_alpha_shape.push_back(pd);
						color_col_alpha_shape.push_back(col);
                    }
			}


			std::vector < Point_d > poin3res;
			for(int k=0; k<(int)poin2res_graham.size();k++){
				FT X1=poin2res_graham[k].x()*vb1.x()+poin2res_graham[k].y()*vb2.x()-plane.d()*vortho.x();
				FT X2=poin2res_graham[k].x()*vb1.y()+poin2res_graham[k].y()*vb2.y()-plane.d()*vortho.y();
				FT X3=poin2res_graham[k].x()*vb1.z()+poin2res_graham[k].y()*vb2.z()-plane.d()*vortho.z();
				Point_d tempo(X1,X2,X3);
				poin3res.push_back(tempo);
			}
			point3dres_graham.push_back(poin3res);
			color_col_graham.push_back(col);

		} 

			char *namo_inliers_colored= "Results_inliers.ply";
			color_inliers(namo_inliers_colored,listpts,vector_of_colors);

			char *namo_alphashape= "Results_alpha_shapes.ply";
			colorplanset2ply(namo_alphashape,point3dres_alpha_shape,color_col_alpha_shape);


return true;
}

#endif
