#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

#include <map>
#include <vector>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;
typedef Polyhedron::Vertex_iterator      Vertex_iterator;
typedef Polyhedron::Vertex_handle        Vertex_handle;
typedef Polyhedron::Edge_iterator        Edge_iterator;
typedef Polyhedron::Halfedge_iterator    Halfedge_iterator;
typedef Polyhedron::Halfedge_handle      Halfedge_handle;
typedef Kernel::Point_3                  Point_3;
typedef Kernel::FT                       FT;
typedef Kernel::Plane_3                  Plane_3;
typedef Kernel::Line_3                   Line_3;
typedef Kernel::Point_2                  Point_2;
typedef Kernel::Vector_3                 Vector_3;

#define PARAM_OFSTREAM_TEXTE std::ios::out
#define OFSTREAM_TEXTE(nomvar, nomfic) std::ofstream nomvar (nomfic, PARAM_OFSTREAM_TEXTE) ; nomvar.setf(std::ios::fixed);

#define RAD2DEG		180.0/M_PI
#define DEG2RAD		M_PI/180.0

#define INF 0xfffffff 

void Dijkstra(const int& s, const int& num_nodes, const std::vector<std::vector<double>>& edge_weights,
      std::vector<int>& flag, std::vector<double>& dis, std::vector<int>& prev);

void convert_edge_to_cylinder(std::string& filename, const std::vector<int>& vertices_on_the_path, std::map<int, Point_3>& map_index_point);

Polyhedron createCylinder(Point_3 p1, Point_3 p2, double size);

Polyhedron createSphere(Point_3 p, double size);

bool save_listpolyhedron2ply(std::vector<Polyhedron> lP, char *nom, std::vector<CGAL::Color> lC, bool with_alpha = false);;

template <class HDS>
class Build_cylinder : public CGAL::Modifier_base<HDS> 
{
private:
	double msize;
	Point_3 mpt1;
	Point_3 mpt2;
	bool mpureTriangle;
	float mstep;
public:
    Build_cylinder(Point_3 pt1,Point_3 pt2, double size, bool pureTriangle, float step = 30.f) {msize = size; mpt1=pt1;  mpt2=pt2; mpureTriangle = pureTriangle; mstep=step;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		const double step = mstep;//works for 30, 45, 90

		int ratiostep = 360/step;
		int nb_point = 2*ratiostep;		//+2 because of the top and bottom points
		int nb_facet = 2*ratiostep;			//nb facet mixed triangle and quad
		//if (mpureTriangle) nb_facet=(nb_facet-2*ratiostep)*2 + 2*ratiostep;

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
       // typedef typename HDS::Vertex   Vertex;

		Line_3 lin(mpt1,mpt2); 
		Plane_3 plane1(mpt1,lin.direction()); 
		Plane_3 plane2(mpt2,lin.direction()); 

		Point_2 ptc=plane1.to_2d(mpt1); 
		std::vector<Point_3> list_pt1;
		std::vector<Point_3> list_pt2;

		for(int i=0;i<ratiostep;i++){
			Point_2 temp_2(cos(i*M_PI*step/180),sin(i*M_PI*step/180)); 
			Vector_3 B1=plane1.base1(); float norme1=sqrt(pow(B1.x(),2)+pow(B1.y(),2)+pow(B1.z(),2)); Vector_3 b1=(1/norme1)*B1;		
			Vector_3 B2=plane1.base2(); float norme2=sqrt(pow(B2.x(),2)+pow(B2.y(),2)+pow(B2.z(),2)); Vector_3 b2=(1/norme2)*B2;	
			Point_3 temp_3(mpt1.x()+msize*(b1.x()*temp_2.x()+b2.x()*temp_2.y()),mpt1.y()+msize*(b1.y()*temp_2.x()+b2.y()*temp_2.y()),mpt1.z()+msize*(b1.z()*temp_2.x()+b2.z()*temp_2.y())); 

			list_pt1.push_back(temp_3);
			list_pt2.push_back(plane2.projection(temp_3));

		}


		//add point
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt1[i]); //std::cout<<list_pt1[i]<<std::endl;}
		for(int i=0;i<ratiostep;i++) B.add_vertex(list_pt2[i]); //std::cout<<list_pt2[i]<<std::endl;}
		
		//add facet 
		for(int i=0;i<ratiostep-1;i++){
			 B.begin_facet();
			 B.add_vertex_to_facet(i);
			 B.add_vertex_to_facet(i+1);
			 B.add_vertex_to_facet(i+ratiostep);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(i+1);
			 B.add_vertex_to_facet(i+ratiostep+1);
			 B.add_vertex_to_facet(i+ratiostep);
			 B.end_facet();
		}
			 
			 B.begin_facet();
			 B.add_vertex_to_facet(ratiostep-1);
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(2*ratiostep-1);
			 B.end_facet();

			 B.begin_facet();
			 B.add_vertex_to_facet(0);
			 B.add_vertex_to_facet(ratiostep);
			 B.add_vertex_to_facet(2*ratiostep-1);
			 B.end_facet();


        B.end_surface();
    }
};

template <class HDS>
class Build_sphere : public CGAL::Modifier_base<HDS> {
private:
	double msize;
	Point_3 mcenter;
	bool mpureTriangle;
public:
    Build_sphere(Point_3 center,double size, bool pureTriangle) {msize = size;mcenter=center;mpureTriangle = pureTriangle;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

		const double step = 30; //works for 30, 45, 90

		int ratiostep = 360/step;
		int nb_point = ratiostep*(ratiostep/2-1)+2;		//+2 because of the top and bottom points
		int nb_facet = ratiostep*(ratiostep/2);			//nb facet mixed triangle and quad
		if (mpureTriangle)
			nb_facet=(nb_facet-2*ratiostep)*2 + 2*ratiostep;

		B.begin_surface( nb_point, nb_facet);			//nb point, nb facet
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( mcenter.x(), mcenter.y(), msize+mcenter.z()));//phi = 90
		double cosphi, sinphi, costheta,sintheta;

		//add point
		double phi;
		int count = 0;
		for (double theta = 0;theta<=360-step;theta+=step)
		{
			
			sintheta = sin(theta*DEG2RAD);costheta = cos(theta*DEG2RAD);
			for (double tp_phi = step;tp_phi<=180-step;tp_phi+=step)//skip phi = 90 and phi = -90
			{
				phi = 90 - tp_phi;count++;
				cosphi = cos(phi*DEG2RAD);sinphi = sin(phi*DEG2RAD);
				B.add_vertex( Point( msize*cosphi*costheta+mcenter.x(), msize*cosphi*sintheta+mcenter.y(), msize*sinphi+mcenter.z()));
			}
		}


		B.add_vertex( Point( mcenter.x(), mcenter.y(), -msize+mcenter.z()));//phi = -90
		
		//add facet (top)
		double step_index = 180/step-1;//= 5
		for (int count = 1;count<=360/step; count++)
		{
        B.begin_facet();
			B.add_vertex_to_facet( 0);
			B.add_vertex_to_facet( (count-1)*step_index+1);

			if (count==360/step)
				B.add_vertex_to_facet(1);
			else
				B.add_vertex_to_facet(count*step_index+1);
        B.end_facet();
		}


		if (mpureTriangle)
		{
			for (int start = 1;start <= step_index-1;start++)
			{
				for (int count = 1;count<=360/step; count++)
				{
				B.begin_facet();
					B.add_vertex_to_facet( start+(count-1)*step_index);
					B.add_vertex_to_facet( start+(count-1)*step_index+1);
					if (count==360/step)
						B.add_vertex_to_facet( start+1 );
					else
						B.add_vertex_to_facet( start+count*step_index+1 );
				B.end_facet();

				B.begin_facet();
					B.add_vertex_to_facet( start+(count-1)*step_index);
					if (count==360/step)
					{
						B.add_vertex_to_facet( start+1 );
						B.add_vertex_to_facet( start);
					}
					else
					{
						B.add_vertex_to_facet( start+count*step_index+1 );
						B.add_vertex_to_facet( start+count*step_index);
					}
				B.end_facet();
				}
			}

		}else{

			for (int start = 1;start <= step_index-1;start++)
			{
				for (int count = 1;count<=360/step; count++)
				{
				B.begin_facet();
					B.add_vertex_to_facet( start+(count-1)*step_index);
					B.add_vertex_to_facet( start+(count-1)*step_index+1);
					if (count==360/step)
					{
						B.add_vertex_to_facet( start+1 );
						B.add_vertex_to_facet( start);
					}
					else
					{
						B.add_vertex_to_facet( start+count*step_index+1 );
						B.add_vertex_to_facet( start+count*step_index);
					}
				B.end_facet();
				}
			}
		}

			
		//add facet (bottom)
		for (int count = 1;count<=360/step; count++)
		{
        B.begin_facet();
			B.add_vertex_to_facet( count*step_index);
			B.add_vertex_to_facet( nb_point-1);

			if (count==360/step)
				B.add_vertex_to_facet(step_index);
			else
				B.add_vertex_to_facet((count+1)*step_index);

        B.end_facet();
		}

	
		

        B.end_surface();
    }
};

int main(int argc, char* argv[])
{
  Polyhedron P;
  std::ifstream in1(argv[1]);
  in1 >> P;

  std::string input_file = argv[1];
  
  std::cout << "number of vertices: " << P.size_of_vertices()
            << ", number of edges: " << P.size_of_halfedges() / 2
            << ", number of facets: " << P.size_of_facets()
            << std::endl;


  std::map<Vertex_handle, int> map_vertex_index;
  std::map<int, Point_3> map_index_point;

  int count_vertex = 0;
  for (Vertex_iterator vi = P.vertices_begin(); vi != P.vertices_end(); vi++) {
    Vertex_handle vh = vi;
    map_index_point.insert(std::make_pair(count_vertex, vh->point()));
    map_vertex_index.insert(std::make_pair(vh, count_vertex++));
  }

  std::cout << "vetex map size: " << map_vertex_index.size() << std::endl;

  // initialize graph edge weight
  int num_nodes = P.size_of_vertices();
  int num_edges = P.size_of_halfedges() / 2;

  std::vector<std::vector<double>> edge_weights;
  for (int i = 0; i < num_nodes; i++) {
    std::vector<double> tmp;
    for (int j = 0; j < num_nodes; j++) {
      if (i == j) {
        tmp.push_back(0.0);
      } else {
        tmp.push_back(std::numeric_limits<double>::max());
      }
    }
    edge_weights.push_back(tmp);
  }

  int count_halfedge = 0;
  for (Halfedge_handle hi = P.halfedges_begin(); hi != P.halfedges_end(); hi++, count_halfedge++) {
    Halfedge_handle h = hi;
    Vertex_handle source = h->vertex();
    Vertex_handle target = h->opposite()->vertex();

    std::map<Vertex_handle, int>::iterator iter = map_vertex_index.find(source);
    if (iter == map_vertex_index.end()) {
      std::cout << "should not be here...";
      return -1;
    }
    int index_source = iter->second;
    Point_3 p1 = iter->first->point();

    iter = map_vertex_index.find(target);
    if (iter == map_vertex_index.end()) {
      std::cout << "should not be here also...";
      return -1;
    }
    int index_target = iter->second;
    Point_3 p2 = iter->first->point();

    // compute edge length
    double distance = std::sqrt(std::pow(p1.x() - p2.x(), 2) + std::pow(p1.y() - p2.y(), 2) + std::pow(p1.z() - p2.z(), 2));

    edge_weights[index_source][index_target] = distance;    
  }

  std::cout << "number of edges: " << count_halfedge / 2 << std::endl;

  // check if edge weight matrix is symmetric
  for (int i = 0; i < num_nodes - 1; i++) {
    for (int j = i + 1; j < num_nodes; j++) {
      // if (edge_weights[i][j] < 100000.0)
      //   std::cout << edge_weights[i][j] << ", ";
      if (edge_weights[i][j] != edge_weights[j][i]) {
        std::cout  << "weight matrix is not symmetric...";
      }
    }
  }

  // iuput index of two vertices
  std::cout << "please input the index of v1 and v2" << std::endl;
  int v1, v2;
  double t1, t2;
  std::cin >> v1 >> v2;
  // while(std::cin >> v1 >> v2) {
    std::cout << "find shortest path from vertex " << v1 << " to " << v2 << std::endl;

    std::vector<int> flag;
    std::vector<double> dis;
    std::vector<int> prev;

    t1 = clock();

    Dijkstra(v1, num_nodes, edge_weights, flag, dis, prev);

    // back propagate
    std::vector<int> vertices_on_the_path;
    vertices_on_the_path.push_back(v2);
    int prev_vertex_index = v2;
    while (prev_vertex_index != v1) {
      prev_vertex_index = prev[prev_vertex_index];
      // std::cout << prev_vertex_index << ", ";
      vertices_on_the_path.push_back(prev_vertex_index);
    }

    std::cout << "Shortest distance from vertex " << v1 << " to " << v2 << " is: " << dis[v2] << std::endl;

		t2 = clock();
		double time = (t2 - t1) / CLOCKS_PER_SEC;
    std::cout << "########" << std::endl;
    std::cout << "Total computing time " << time << " sec." << std::endl;
    std::cout << "########" << std::endl;

    // check
    double dist = 0.0;
    for (int i = 0; i < vertices_on_the_path.size() - 1; i++) {
      // std::cout << i << "," << vertices_on_the_path[i] << ", " << vertices_on_the_path[i+1] << std::endl;
      Point_3 p1 = map_index_point[vertices_on_the_path[i]];
      Point_3 p2 = map_index_point[vertices_on_the_path[i+1]];
      double distance = std::sqrt(std::pow(p1.x() - p2.x(), 2) + std::pow(p1.y() - p2.y(), 2) + std::pow(p1.z() - p2.z(), 2));
      dist += distance;
    }

    std::cout << "check distance: " << dist << std::endl;

    // output result
    convert_edge_to_cylinder(input_file, vertices_on_the_path, map_index_point);
  // }

  return -1;
}

void Dijkstra(const int& s, const int& num_nodes, const std::vector<std::vector<double>>& edge_weights,
    std::vector<int>& flag, std::vector<double>& dis, std::vector<int>& prev) {
  	// TO START FROM HERE (YOU)
	//-----------------------
	int i,j,k;
	//初始化所有节点
	flag[0] = {0};
	//第一个节点已经访问
	flag[s] = 1;
	//初始化源点到其他点的距离
	for(i = 1; i<= num_nodes; i++)
	{
		dis[i] = edge_weights[s][i];
		if(dis[i] == INF || dis[i] == 0)
		{
			prev[i] = 0;
		}
		else
		{
			prev[i] = s;
		}
	}
	int min;
	for (i = 1; i<= num_nodes; i++)
	{
		min = INF;
		for (j=1; j <= num_nodes; j++)
		{
			if(!flag[j] && dis[j] < min)
			{
				min = dis[j];
				k = j;
			}
		}

		flag[k] = 1;

		for (j = 1; j <= num_nodes; j++)
		{
			if (dis[j] > dis[k] + edge_weights[k][i])
			{
				dis[j] = dis[k] + edge_weights[k][j];
				prev[j] = k;
			}
		}
	}

	//-----------------------
  return;
}

void convert_edge_to_cylinder(std::string& filename, const std::vector<int>& vertices_on_the_path, std::map<int, Point_3>& map_index_point) {
	std::vector<Point_3> vertices;
	
	// 2. for each edge, compute a cyclinder
	std::vector<Polyhedron> P;
	std::vector<CGAL::Color> C;
	double radius_cylinder = 0.001;
	
	for (int i = 0; i < vertices_on_the_path.size() - 1; ++i) {
		Point_3 pt1 = map_index_point[vertices_on_the_path[i]];
		Point_3 pt2 = map_index_point[vertices_on_the_path[i+1]];

		Polyhedron P_cylinder = createCylinder(pt1, pt2, radius_cylinder);
		P.push_back(P_cylinder);

		//compute color

		//int r = 255, g = 255, b = 102; // yellow
		int r = 255, g = 102, b = 102; // red
		//int r = 102, g = 178, b = 255; // blue
		//int r = 32, g = 32, b = 32; // black

		CGAL::Color c(r, g, b);
		C.push_back(c);
	}

	Point_3 pt1 = map_index_point[vertices_on_the_path[0]];
	Polyhedron P_sphere_1 = createSphere(pt1, 5 * radius_cylinder);

	P.push_back(P_sphere_1);

	//compute color

	int r = 255, g = 255, b = 102; // yellow
	//int r = 255, g = 102, b = 102; // red
	//int r = 102, g = 178, b = 255; // blue
	//int r = 32, g = 32, b = 32; // black

	CGAL::Color c(r, g, b);
	C.push_back(c);

	Point_3 pt2 = map_index_point[vertices_on_the_path.back()];
	Polyhedron P_sphere_2 = createSphere(pt2, 5 * radius_cylinder);

	P.push_back(P_sphere_2);
	C.push_back(c);

	// 3. output cylinders to PLY file
	std::string output_filename = filename.substr(0, filename.find_last_of('.')) + "_edges_on_the_shortest_path.ply";
	char* out_filename = new char[output_filename.length() + 1];
	strcpy(out_filename, output_filename.c_str());

	save_listpolyhedron2ply(P, out_filename, C);
}

Polyhedron createCylinder(Point_3 p1, Point_3 p2, double size) {
	Polyhedron P;
	Build_cylinder<Polyhedron::HalfedgeDS> poly1(p1,p2,size,true);
	P.delegate(poly1);

	return P;
}

Polyhedron createSphere(Point_3 p, double size) {
	Polyhedron P;
	Build_sphere<Polyhedron::HalfedgeDS> poly1(p,size,true);
	P.delegate(poly1);
	return P;
}

bool save_listpolyhedron2ply(std::vector<Polyhedron> lP, char *nom, std::vector<CGAL::Color> lC, bool with_alpha) {
	std::cout <<"Saving polyhedron to "<<nom<<std::endl;
	OFSTREAM_TEXTE(fic,nom);
	
	int nb_vertices = 0;
	int nb_facets = 0;

	std::vector<Polyhedron>::iterator it;
	for(it = lP.begin(); it != lP.end(); it++)
	{				
		nb_vertices+=(*it).size_of_vertices();
		nb_facets+=(*it).size_of_facets ();
	}

	//header
	fic << "ply" << std::endl;
	fic << "format ascii 1.0" << std::endl;
	fic << "element vertex " << nb_vertices << std::endl;
	fic << "property float x" << std::endl;
	fic << "property float y" << std::endl;
	fic << "property float z" << std::endl;
	fic << "element face " << nb_facets << std::endl;
	fic << "property list uchar int vertex_index" << std::endl;
	fic << "property uchar red" << std::endl;
	fic << "property uchar green" << std::endl;
	fic << "property uchar blue" << std::endl;
	if (with_alpha)
		fic << "property uchar alpha" << std::endl;
	fic << "end_header" << std::endl;

	//vertex list
	std::map<Point_3, int> map_point_index;
	int count = 0;
	for (it = lP.begin(); it != lP.end(); it++)
	{
		for (Polyhedron::Vertex_const_iterator vi = (*it).vertices_begin(); vi != (*it).vertices_end(); ++vi)
		{
			/*std::map<Point_d, int>::iterator iter = map_point_index.find(vi->point());
			if (iter == map_point_index.end())
			{*/
			map_point_index[vi->point()] = count++;
			fic << vi->point().x() << " " << vi->point().y() << " " << vi->point().z() << std::endl;
			//}
		}
	}
	

	int i=0;
	for(it = lP.begin(); it != lP.end(); it++, i++)
	{	
		//facet list
		Polyhedron::Halfedge_around_facet_const_circulator hc, hc_end;
		for(Polyhedron::Facet_const_iterator fi = (*it).facets_begin(); fi != (*it).facets_end(); ++fi)
		{
			hc = fi->facet_begin();
			hc_end = hc;
			std::size_t n = circulator_size( hc);
			CGAL_assertion( n >= 3);
			fic<<n<<" ";
			do {
				fic<<map_point_index[hc->vertex()->point()]<<" ";
				++hc;
			} while( hc != hc_end);
			fic<<(int)lC[i].r()<<" "<<(int)lC[i].g()<<" "<<(int)lC[i].b();
			if (with_alpha) fic<<" "<<(int)lC[i].alpha();
			fic<<std::endl;
		}
	}


		
fic<<std::endl;	

return true;
}

