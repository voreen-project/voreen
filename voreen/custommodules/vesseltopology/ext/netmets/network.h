#ifndef STIM_NETWORK_H
#define STIM_NETWORK_H

#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <stim/math/vec3.h>
#include <stim/visualization/obj.h>
#include <stim/visualization/swc.h>
#include "cylinder.h"
#include <stim/cuda/cudatools/timer.h>
#include <stim/cuda/cudatools/callable.h>
#include <stim/structures/kdtree.cuh>
//********************help function********************
// gaussian_function
CUDA_CALLABLE float gaussianFunction(float x, float std = 25) { return exp(-x / (2 * std*std)); }  // std default sigma value is 25

// compute metric in parallel
#ifdef __CUDACC__
template <typename T>
__global__ void find_metric_parallel(T* M, size_t n, T* D, float sigma){
	size_t x = blockDim.x * blockIdx.x + threadIdx.x;
	if(x >= n) return;
	M[x] = 1.0f - gaussianFunction(D[x], sigma);
}

//find the corresponding edge index from array index
__global__ void find_edge_index_parallel(size_t* I, size_t n, unsigned* R, size_t* E, size_t ne){
	size_t x = blockDim.x * blockIdx.x + threadIdx.x;
	if(x >= n) return;
	unsigned i = 0;
	size_t N = 0;
	for(unsigned e = 0; e < ne; e++){
		N += E[e];
		if(I[x] < N){
			R[x] = i;
			break;
		}
		i++;
	}
}
#endif

//hard-coded factor
int threshold_fac;

namespace stim{
/** This is the a class that interfaces with gl_spider in order to store the currently
 *   segmented network. The following data is stored and can be extracted:
 *   1)Network geometry and centerline.
 *   2)Network connectivity (a graph of nodes and edges), reconstructed using kdtree.
*/

template<typename T>
class network{

	///Each edge is a fiber with two nodes.
	///Each node is an in index to the endpoint of the fiber in the nodes array.
	class edge : public cylinder<T>
	{
		public:
	
		unsigned int v[2];		//unique id's designating the starting and ending
		// default constructor
		edge() : cylinder<T>() {
			v[1] = (unsigned)(-1); v[0] = (unsigned)(-1);
		}
		/// Constructor - creates an edge from a list of points by calling the stim::fiber constructor
/*
		///@param v0, the starting index.
		///@param v1, the ending index.
		///@param sz, the number of point in the fiber.
		edge(unsigned int v0, unsigned int v1, unsigned int sz) : cylinder<T>(
		{

		}
*/
		edge(std::vector<stim::vec3<T> > p, std::vector<T> s)
			: cylinder<T>(p,s)
		{
		}
		///@param p is an array of positions in space
		edge(stim::centerline<T> p) : cylinder<T>(p){}

		/// Copy constructor creates an edge from a cylinder
		edge(stim::cylinder<T> f) : cylinder<T>(f) {}

		/// Resamples an edge by calling the fiber resampling function
		edge resample(T spacing){
			edge e(cylinder<T>::resample(spacing));	//call the fiber->edge constructor
			e.v[0] = v[0];					//copy the vertex data
			e.v[1] = v[1];

			return e;						//return the new edge
		}

		/// Output the edge information as a string
		std::string str(){
			std::stringstream ss;
			ss<<"("<<cylinder<T>::size()<<")\tl = "<<this->length()<<"\t"<<v[0]<<"----"<<v[1];
			return ss.str();
		}

		std::vector<edge> split(unsigned int idx){
		
			std::vector< stim::cylinder<T> > C;
			C.resize(2);
			C =	(*this).cylinder<T>::split(idx);
			std::vector<edge> E(C.size());

			for(unsigned e = 0; e < E.size(); e++){
				E[e] = C[e];
			}
			return E;
		}

		/// operator for writing the edge information into a binary .nwt file.
		friend std::ofstream& operator<<(std::ofstream& out, const edge& e)
		{
			out.write(reinterpret_cast<const char*>(&e.v[0]), sizeof(unsigned int));	///write the starting point.
			out.write(reinterpret_cast<const char*>(&e.v[1]), sizeof(unsigned int));	///write the ending point.
			unsigned int sz = e.size();	///write the number of point in the edge.
			out.write(reinterpret_cast<const char*>(&sz), sizeof(unsigned int));
			for(int i = 0; i < sz; i++)	///write each point
			{
				stim::vec3<T> point = e[i];
				out.write(reinterpret_cast<const char*>(&point[0]), 3*sizeof(T));
			//	for(int j = 0; j < nmags(); j++)	//future code for multiple mags
			//	{
				out.write(reinterpret_cast<const char*>(&e.R[i]), sizeof(T));	///write the radius
				//std::cout << point.str() << " " << e.R[i] << std::endl;
			//	}
			}
			return out;	//return stream
		}

		/// operator for reading an edge from a binary .nwt file.
		friend std::ifstream& operator>>(std::ifstream& in, edge& e)
		{
			unsigned int v0, v1, sz;
			in.read(reinterpret_cast<char*>(&v0), sizeof(unsigned int));	//read the staring point.
			in.read(reinterpret_cast<char*>(&v1), sizeof(unsigned int));	//read the ending point
			in.read(reinterpret_cast<char*>(&sz), sizeof(unsigned int));	//read the number of points in the edge
//			stim::centerline<T> temp = stim::centerline<T>(sz);		//allocate the new edge
//			e = edge(temp);
			std::vector<stim::vec3<T> > p(sz);
			std::vector<T> r(sz);
			for(int i = 0; i < sz; i++)		//set the points and radii to the newly read values
			{
				stim::vec3<T> point;
				in.read(reinterpret_cast<char*>(&point[0]), 3*sizeof(T));
				p[i] = point;
				T mag;
		//				for(int j = 0; j < nmags(); j++)		///future code for mags
		//				{
				in.read(reinterpret_cast<char*>(&mag), sizeof(T));	
				r[i] = mag;
				//std::cout << point.str() << " " << mag << std::endl;
		//				}
			}
			e = edge(p,r);
			e.v[0] = v0; e.v[1] = v1;
			return in;
		}
	};

	///Node class that stores the physical position of the node as well as the edges it is connected to (edges that connect to it), As well as any additional data necessary.
	class vertex : public stim::vec3<T>
	{
		public:
			//std::vector<unsigned int> edges;					//indices of edges connected to this node.
			std::vector<unsigned int> e[2];						//indices of edges going out (e[0]) and coming in (e[1])
			//stim::vec3<T> p;							//position of this node in physical space.
			//default constructor
			vertex() : stim::vec3<T>()
			{
			}
			//constructor takes a stim::vec
			vertex(stim::vec3<T> p) : stim::vec3<T>(p){}

			/// Output the vertex information as a string
			std::string 
			str(){
				std::stringstream ss;
				ss<<"\t(x, y, z) = "<<stim::vec3<T>::str();

				if(e[0].size() > 0){
					ss<<"\t> ";
					for(unsigned int o = 0; o < e[0].size(); o++)
						ss<<e[0][o]<<" ";
				}
				if(e[1].size() > 0){
					ss<<"\t< ";
					for(unsigned int i = 0; i < e[1].size(); i++)
						ss<<e[1][i]<<" ";
				}

				return ss.str();
			}
			///operator for writing the vector into the stream;
			friend std::ofstream& operator<<(std::ofstream& out, const vertex& v)
			{
				unsigned int s0, s1;
				s0 = v.e[0].size();
				s1 = v.e[1].size();
				out.write(reinterpret_cast<const char*>(&v.ptr[0]), 3*sizeof(T));	///write physical vertex location
				out.write(reinterpret_cast<const char*>(&s0), sizeof(unsigned int));	///write the number of "outgoing edges"
				out.write(reinterpret_cast<const char*>(&s1), sizeof(unsigned int));	///write the number of "incoming edges"	
				if (s0 != 0)
					out.write(reinterpret_cast<const char*>(&v.e[0][0]), sizeof(unsigned int)*v.e[0].size());	///write the "outgoing edges"
				if (s1 != 0)
					out.write(reinterpret_cast<const char*>(&v.e[1][0]), sizeof(unsigned int)*v.e[1].size());	///write the "incoming edges"
				return out;
			}

			///operator for reading the vector out of the stream;
			friend std::ifstream& operator>>(std::ifstream& in, vertex& v)
			{
				in.read(reinterpret_cast<char*>(&v[0]), 3*sizeof(T));	///read the physical position
				unsigned int s[2];					
				in.read(reinterpret_cast<char*>(&s[0]), 2*sizeof(unsigned int));	///read the sizes of incoming and outgoing edge arrays

				std::vector<unsigned int> one(s[0]);
				std::vector<unsigned int> two(s[1]);
				v.e[0] = one;
				v.e[1] = two;
				if (one.size() != 0)
					in.read(reinterpret_cast<char*>(&v.e[0][0]), s[0] * sizeof(unsigned int));		///read the arrays of "outgoing edges"
				if (two.size() != 0)
					in.read(reinterpret_cast<char*>(&v.e[1][0]), s[1] * sizeof(unsigned int));		///read the arrays of "incoming edges"
				return in;
			}

	};

protected:

	std::vector<edge> E;       //list of edges
	std::vector<vertex> V;	   //list of vertices.

public:

	///default constructor
	network()
	{
		
	}

	///constructor with a file to load.
	network(std::string fileLocation)
	{
		load_obj(fileLocation);
	}

	///Returns the number of edges in the network.
	unsigned int edges(){
		return E.size();
	}

	///Returns the number of nodes in the network.
	unsigned int vertices(){
		return V.size();
	}

	///Returns the radius at specific point in the edge
	T get_r(unsigned e, unsigned i) {
		return E[e].r(i);
	}

	///Returns the average radius of specific edge
	T get_average_r(unsigned e) {
		T result = 0.0;
		unsigned n = E[e].size();
		for (unsigned p = 0; p < n; p++)
			result += E[e].r(p);

		return (T)result / n;
	}

	///Returns the length of current edge
	T get_l(unsigned e) {
		return E[e].length();
	}

	///Returns the start vertex of current edge
	size_t get_start_vertex(unsigned e) {
		return E[e].v[0];
	}

	///Returns the end vertex of current edge
	size_t get_end_vertex(unsigned e) {
		return E[e].v[1];
	}

	///Returns one vertex
	stim::vec3<T> get_vertex(unsigned i) {
		return V[i];
	}

	///Returns the boundary vertices' indices
	std::vector<unsigned> get_boundary_vertex() {
		std::vector<unsigned> result;

		for (unsigned v = 0; v < V.size(); v++) {
			if (V[v].e[0].size() + V[v].e[1].size() == 1) {	// boundary vertex
				result.push_back(v);
			}
		}

		return result;
	}

	///Set radius
	void set_r(unsigned e, std::vector<T> radius) {
		E[e].cylinder<T>::copy_r(radius);
	}

	void set_r(unsigned e, T radius) {
		for (size_t i = 0; i < E[e].size(); i++)
			E[e].cylinder<T>::set_r(i, radius);
	}
	//scale the network by some constant value
	//	I don't think these work??????
	/*std::vector<vertex> operator*(T s){
		for (unsigned i=0; i< vertices; i ++ ){
			V[i] = V[i] * s;
		}
		return V;
	}

	std::vector<vertex> operator*(vec<T> s){
		for (unsigned i=0; i< vertices; i ++ ){
			for (unsigned dim = 0 ; dim< 3; dim ++){
				V[i][dim] = V[i][dim] * s[dim];
			}
		}
		return V;
	}*/

	// Returns an average of branching index in the network
	double BranchingIndex(){
		double B=0;
		for(unsigned v=0; v < V.size(); v ++){
			B += ((V[v].e[0].size()) + (V[v].e[1].size()));
		}
		B = B / V.size();
		return B;

	}

	// Returns number of branch points in thenetwork
	unsigned int BranchP(){
		unsigned int B=0;
		unsigned int c;
		for(unsigned v=0; v < V.size(); v ++){
			c = ((V[v].e[0].size()) + (V[v].e[1].size()));
			if (c > 2){
			B += 1;}
		}		
		return B;

	}

	// Returns number of end points (tips) in thenetwork
	unsigned int EndP(){
		unsigned int B=0;
		unsigned int c;
		for(unsigned v=0; v < V.size(); v ++){
			c = ((V[v].e[0].size()) + (V[v].e[1].size()));
			if (c == 1){
			B += 1;}
		}		
		return B;

	}

	//// Returns a dictionary with the key as the vertex
	//std::map<std::vector<vertex>,unsigned int> DegreeDict(){
	//	std::map<std::vector<vertex>,unsigned int> dd;
	//	unsigned int c = 0;
	//	for(unsigned v=0; v < V.size(); v ++){
	//		c = ((V[v].e[0].size()) + (V[v].e[1].size()));
	//		dd[V[v]] = c;
	//	}
	//	return dd;
	//}

	//// Return number of branching stems
	//unsigned int Stems(){
	//	unsigned int s = 0;
	//	std::map<std::vector<vertex>,unsigned int> dd;
	//	dd = DegreeDict();
	//	//for(unsigned v=0; v < V.size(); v ++){
	//	//	V[v].e[0].
	//	return s;
	//}

	//Calculate Metrics---------------------------------------------------
	// Returns an average of fiber/edge lengths in the network
	double Lengths(){
		stim::vec<T> L;
		double sumLength = 0;
		for(unsigned e = 0; e < E.size(); e++){				//for each edge in the network
			L.push_back(E[e].length());						//append the edge length
			sumLength = sumLength + E[e].length();
		}
		double avg = sumLength / E.size();
		return avg;
	}


	// Returns an average of tortuosities in the network
	double Tortuosities(){
		stim::vec<T> t;
		stim::vec<T> id1, id2;                        // starting and ending vertices of the edge
		double distance;double tortuosity;double sumTortuosity = 0;
		for(unsigned e = 0; e < E.size(); e++){				//for each edge in the network
			id1 = E[e][0];									//get the edge starting point
			id2 = E[e][E[e].size() - 1];					//get the edge ending point
			distance = (id1 - id2).len();                   //displacement between the starting and ending points
			if(distance > 0){
				tortuosity = E[e].length()/	distance	;		// tortuoisty = edge length / edge displacement
			}
			else{
				tortuosity = 0;}
			t.push_back(tortuosity);
			sumTortuosity += tortuosity;
		}
		double avg = sumTortuosity / E.size();
		return avg;
	}

	// Returns average contraction of the network
	double Contractions(){
	stim::vec<T> t;
	stim::vec<T> id1, id2;                        // starting and ending vertices of the edge
	double distance;double contraction;double sumContraction = 0;
	for(unsigned e = 0; e < E.size(); e++){				//for each edge in the network
		id1 = E[e][0];									//get the edge starting point
		id2 = E[e][E[e].size() - 1];					//get the edge ending point
		distance = (id1 - id2).len();                   //displacement between the starting and ending points
		contraction = distance / E[e].length();		// tortuoisty = edge length / edge displacement
		t.push_back(contraction);
		sumContraction += contraction;
	}
	double avg = sumContraction / E.size();
	return avg;
	}

	// returns average fractal dimension of the branches of the network
	double FractalDimensions(){
	stim::vec<T> t;
	stim::vec<T> id1, id2;                        // starting and ending vertices of the edge
	double distance;double fract;double sumFractDim = 0;
	for(unsigned e = 0; e < E.size(); e++){				//for each edge in the network
		id1 = E[e][0];									//get the edge starting point
		id2 = E[e][E[e].size() - 1];					//get the edge ending point
		distance = (id1 - id2).len();                   //displacement between the starting and ending points
		fract = std::log(distance) / std::log(E[e].length());		// tortuoisty = edge length / edge displacement
		t.push_back(sumFractDim);
		sumFractDim += fract;
	}
	double avg = sumFractDim / E.size();
	return avg;
	}

	//returns a cylinder represented a given fiber (based on edge index)
	stim::cylinder<T> get_cylinder(unsigned e){
		return E[e];									//return the specified edge (casting it to a fiber)
	}

	/// subdivide current network
	void subdivision() {

		std::vector<unsigned> ori_index;		// original index
		std::vector<unsigned> new_index;		// new index
		std::vector<edge> nE;					// new edge
		std::vector<vertex> nV;					// new vector
		unsigned id = 0;
		unsigned num_edge = (*this).E.size();

		for (unsigned i = 0; i < num_edge; i++) {
			if (E[i].size() == 2) {				// if current edge can't be subdivided
				stim::centerline<T> line(2);
				for (unsigned k = 0; k < 2; k++)
					line[k] = E[i][k];
				line.update();

				edge new_edge(line);

				vertex new_vertex = new_edge[0];
				id = E[i].v[0];
				auto position = std::find(ori_index.begin(), ori_index.end(), id);
				if (position == ori_index.end()) {		 // new vertex
					ori_index.push_back(id);
					new_index.push_back(nV.size());

					new_vertex.e[0].push_back(nE.size());
					new_edge.v[0] = nV.size();
					nV.push_back(new_vertex);			// push back vertex as a new vertex
				}
				else {									// existing vertex
					int k = std::distance(ori_index.begin(), position);
					new_edge.v[0] = new_index[k];
					nV[new_index[k]].e[0].push_back(nE.size());
				}

				new_vertex = new_edge[1];
				id = E[i].v[1];
				position = std::find(ori_index.begin(), ori_index.end(), id);
				if (position == ori_index.end()) {		 // new vertex
					ori_index.push_back(id);
					new_index.push_back(nV.size());

					new_vertex.e[1].push_back(nE.size());
					new_edge.v[1] = nV.size();
					nV.push_back(new_vertex);			// push back vertex as a new vertex
				}
				else {									// existing vertex
					int k = std::distance(ori_index.begin(), position);
					new_edge.v[1] = new_index[k];
					nV[new_index[k]].e[1].push_back(nE.size());
				}

				nE.push_back(new_edge);

				nE[nE.size() - 1].cylinder<T>::set_r(0, E[i].cylinder<T>::r(0));
				nE[nE.size() - 1].cylinder<T>::set_r(1, E[i].cylinder<T>::r(1));
			}
			else {								// subdivide current edge
				for (unsigned j = 0; j < E[i].size() - 1; j++) {
					stim::centerline<T> line(2);
					for (unsigned k = 0; k < 2; k++)
						line[k] = E[i][j + k];
					line.update();

					edge new_edge(line);

					if (j == 0) {						// edge contains original starting point
						vertex new_vertex = new_edge[0];
						id = E[i].v[0];
						auto position = std::find(ori_index.begin(), ori_index.end(), id);
						if (position == ori_index.end()) {		 // new vertex
							ori_index.push_back(id);
							new_index.push_back(nV.size());

							new_vertex.e[0].push_back(nE.size());
							new_edge.v[0] = nV.size();
							nV.push_back(new_vertex);			// push back vertex as a new vertex
						}
						else {									// existing vertex
							int k = std::distance(ori_index.begin(), position);
							new_edge.v[0] = new_index[k];
							nV[new_index[k]].e[0].push_back(nE.size());
						}

						new_vertex = new_edge[1];
						new_vertex.e[1].push_back(nE.size());
						new_edge.v[1] = nV.size();
						nV.push_back(new_vertex);				// push back internal point as a new vertex

						nE.push_back(new_edge);
					}

					else if (j == E[i].size() - 2) {	// edge contains original ending point

						vertex new_vertex = new_edge[1];
						nV[nV.size() - 1].e[0].push_back(nE.size());
						new_edge.v[0] = nV.size() - 1;

						id = E[i].v[1];
						auto position = std::find(ori_index.begin(), ori_index.end(), id);
						if (position == ori_index.end()) {		 // new vertex
							ori_index.push_back(id);
							new_index.push_back(nV.size());

							new_vertex.e[1].push_back(nE.size());
							new_edge.v[1] = nV.size();
							nV.push_back(new_vertex);			// push back vertex as a new vertex
						}
						else {									// existing vertex
							int k = std::distance(ori_index.begin(), position);
							new_edge.v[1] = new_index[k];
							nV[new_index[k]].e[1].push_back(nE.size());
						}

						nE.push_back(new_edge);
					}

					else {
						vertex new_vertex = new_edge[1];

						nV[nV.size() - 1].e[0].push_back(nE.size());
						new_vertex.e[1].push_back(nE.size());
						new_edge.v[0] = nV.size() - 1;
						new_edge.v[1] = nV.size();
						nV.push_back(new_vertex);

						nE.push_back(new_edge);
					}

					// get radii
					nE[nE.size() - 1].cylinder<T>::set_r(0, E[i].cylinder<T>::r(j));
					nE[nE.size() - 1].cylinder<T>::set_r(1, E[i].cylinder<T>::r(j + 1));
				}
			}
		}

		(*this).E = nE;
		(*this).V = nV;
	}

	//load a network from an OBJ file
	void load_obj(std::string filename){

		stim::obj<T> O;									//create an OBJ object
		O.load(filename);								//load the OBJ file as an object

		std::vector<unsigned> id2vert;							//this list stores the OBJ vertex ID associated with each network vertex

		unsigned i[2];									//temporary, IDs associated with the first and last points in an OBJ line

		//for each line in the OBJ object
		for(unsigned int l = 1; l <= O.numL(); l++){

			std::vector< stim::vec<T> > c;						//allocate an array of points for the vessel centerline
			O.getLine(l, c);							//get the fiber centerline

			stim::centerline<T> c3(c.size());
			for(size_t j = 0; j < c.size(); j++)
				c3[j] = c[j];
			c3.update();

	//		edge new_edge = c3;		///This is dangerous.
			edge new_edge(c3);
					
			//create an edge from the given centerline
			unsigned int I = new_edge.size();					//calculate the number of points on the centerline

			//get the first and last vertex IDs for the line
			std::vector< unsigned > id;						//create an array to store the centerline point IDs
			O.getLinei(l, id);							//get the list of point IDs for the line
			i[0] = id.front();							//get the OBJ ID for the first element of the line
			i[1] = id.back();							//get the OBJ ID for the last element of the line

			std::vector<unsigned>::iterator it;					//create an iterator for searching the id2vert array
			unsigned it_idx;							//create an integer for the id2vert entry index

			//find out if the nodes for this fiber have already been created
			it = find(id2vert.begin(), id2vert.end(), i[0]);	//look for the first node
			if(it == id2vert.end()){							//if i[0] hasn't already been used
				vertex new_vertex = new_edge[0];				//create a new vertex, assign it a position
				bool flag = false;
				unsigned j = 0;
				for (; j < V.size(); j++) {						// check whether current vertex is already exist
					if (new_vertex == V[j]) {
						flag = true;
						break;
					}
				}
				if (!flag) {									// unique one
					new_vertex.e[0].push_back(E.size());				//add the current edge as outgoing
					new_edge.v[0] = V.size();					//add the new edge to the edge
					V.push_back(new_vertex);					//add the new vertex to the vertex list
					id2vert.push_back(i[0]);					//add the ID to the ID->vertex conversion list
				}
				else {
					V[j].e[0].push_back(E.size());
					new_edge.v[0] = j;
				}
			}
			else{									//if the vertex already exists
				it_idx = std::distance(id2vert.begin(), it);
				V[it_idx].e[0].push_back(E.size());				//add the current edge as outgoing
				new_edge.v[0] = it_idx;
			}

			it = find(id2vert.begin(), id2vert.end(), i[1]);			//look for the second ID
			if(it == id2vert.end()){						//if i[1] hasn't already been used
				vertex new_vertex = new_edge[I-1];				//create a new vertex, assign it a position
				bool flag = false;
				unsigned j = 0;
				for (; j < V.size(); j++) {					// check whether current vertex is already exist
					if (new_vertex == V[j]) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					new_vertex.e[1].push_back(E.size());				//add the current edge as incoming
					new_edge.v[1] = V.size();                                  	//add the new vertex to the edge
					V.push_back(new_vertex);					//add the new vertex to the vertex list
					id2vert.push_back(i[1]);					//add the ID to the ID->vertex conversion list
				}
				else {
					V[j].e[1].push_back(E.size());
					new_edge.v[1] = j;
				}
			}
			else{									//if the vertex already exists
				it_idx = std::distance(id2vert.begin(), it);
				V[it_idx].e[1].push_back(E.size());				//add the current edge as incoming
				new_edge.v[1] = it_idx;
			}

			E.push_back(new_edge);							//push the edge to the list

		}

		// copy the radii information from OBJ
		/*if (O.numVT()) {
			unsigned k = 0;
			for (unsigned i = 0; i < E.size(); i++) {
				for (unsigned j = 0; j < E[i].size(); j++) {
					E[i].cylinder<T>::set_r(j, O.getVT(k)[0] / 2);
					k++;
				}
			}
		}*/
		// OBJ class assumes that in L the two values are equal
		if (O.numVT()) {
			std::vector< unsigned > id;						//create an array to store the centerline point IDs
			for (unsigned i = 0; i < O.numL(); i++) {
				id.clear();
				O.getLinei(i + 1, id);							//get the list of point IDs for the line
				for (unsigned j = 0; j < id.size(); j++)
					E[i].cylinder<T>::set_r(j, O.getVT(id[j] - 1)[0] / 2);
			}
		}
	}

	///loads a .nwt file. Reads the header and loads the data into the network according to the header.
	void
	loadNwt(std::string filename)
	{
		int dims[2];		///number of vertex, number of edges
		readHeader(filename, &dims[0]);		//read header
		std::ifstream file;
		file.open(filename.c_str(), std::ios::in | std::ios::binary);		///skip header information.
		file.seekg(14+58+4+4, file.beg);
		vertex v;
		for(int i = 0; i < dims[0]; i++)		///for every vertex, read vertex, add to network.
		{
			file >> v;
			V.push_back(v);
//			std::cout << i << " " << v.str() << std::endl;
		}

		std::cout << std::endl;
		for(int i = 0; i < dims[1]; i++)		///for every edge, read edge, add to network.
		{
			edge e;
			file >> e;
			E.push_back(e);
			//std::cout << i << " " << E[i].str() << std::endl;		// not necessary?
		}
		file.close();
	}

	///saves a .nwt file. Writes the header in raw text format, then saves the network as a binary file.
	void
	saveNwt(std::string filename)
	{
		writeHeader(filename);
		std::ofstream file;
		file.open(filename.c_str(), std::ios::out | std::ios::binary | std::ios::app);	///since we have written the header we are not appending.
		for(int i = 0; i < V.size(); i++)	///look through the Vertices and write each one.
		{
//			std::cout << i << " " << V[i].str() << std::endl;
			file << V[i];
		}
		for(int i = 0; i < E.size(); i++)	///loop through the Edges and write each one.
		{
			//std::cout << i << " " << E[i].str() << std::endl;		// not necesarry?
			file << E[i];
		}
		file.close();
	}


	///Writes the header information to a .nwt file.
	void
	writeHeader(std::string filename)
	{
		std::string magicString = "nwtFileFormat ";		///identifier for the file.
		std::string desc = "fileid(14B), desc(58B), #vertices(4B), #edges(4B): bindata";
		int hNumVertices = V.size();		///int byte header storing the number of vertices in the file
		int hNumEdges = E.size();		///int byte header storing the number of edges.
		std::ofstream file;
		file.open(filename.c_str(), std::ios::out | std::ios::binary);
		std::cout << hNumVertices << " " << hNumEdges << std::endl;
		file.write(reinterpret_cast<const char*>(&magicString.c_str()[0]), 14);	//write the file id
		file.write(reinterpret_cast<const char*>(&desc.c_str()[0]), 58);	//write the description
		file.write(reinterpret_cast<const char*>(&hNumVertices), sizeof(int));	//write #vert.
		file.write(reinterpret_cast<const char*>(&hNumEdges), sizeof(int));	//write #edges
//		file << magicString.c_str() << desc.c_str() << hNumVertices << hNumEdges;
		file.close();
		
	}

	///Reads the header information from a .nwt file.
	void
	readHeader(std::string filename, int *dims)
	{
		char magicString[14];		///id
		char desc[58];			///description
		int hNumVertices;		///#vert
		int hNumEdges;			///#edges
		std::ifstream file;		////create stream
		file.open(filename.c_str(), std::ios::in | std::ios::binary);
		file.read(reinterpret_cast<char*>(&magicString[0]), 14);	///read the file id.
		file.read(reinterpret_cast<char*>(&desc[0]), 58);		///read the description
		file.read(reinterpret_cast<char*>(&hNumVertices), sizeof(int));	///read the number of vertices
		file.read(reinterpret_cast<char*>(&hNumEdges), sizeof(int));	///read the number of edges
//		std::cout << magicString << desc << hNumVertices << " " <<  hNumEdges << std::endl;
		file.close();							///close the file.
		dims[0] = hNumVertices;						///fill the returned reference.
		dims[1] = hNumEdges;
	}

	//load a network from an SWC file
	void load_swc(std::string filename) {
		stim::swc<T> S;										// create swc variable
		S.load(filename);									// load the node information
		S.create_tree();									// link those node according to their linking relationships as a tree
		S.resample();

		//NT.push_back(S.node[0].type);						// set the neuronal_type value to the first vertex in the network
		std::vector<unsigned> id2vert;						// this list stores the SWC vertex ID associated with each network vertex
		unsigned i[2];										// temporary, IDs associated with the first and last points

		for (unsigned int l = 0; l < S.numE(); l++) {		// for every edge
			//NT.push_back(S.node[l].type);

			std::vector< stim::vec3<T> > c;
			S.get_points(l, c);

			stim::centerline<T> c3(c.size());				// new fiber
			
			for (unsigned j = 0; j < c.size(); j++)
				c3[j] = c[j];								// copy the points
		
			c3.update();									// update the L information
			
			stim::cylinder<T> C3(c3);						// create a new cylinder in order to copy the origin radius information
			// upadate the R information
			std::vector<T> radius;
			S.get_radius(l, radius);

			C3.copy_r(radius);

			edge new_edge(C3);								// new edge	

			//create an edge from the given centerline
			unsigned int I = (unsigned)new_edge.size();				//calculate the number of points on the centerline
			
			//get the first and last vertex IDs for the line
			i[0] = S.E[l].front();
			i[1] = S.E[l].back();

			std::vector<unsigned>::iterator it;				//create an iterator for searching the id2vert array
			unsigned it_idx;								//create an integer for the id2vert entry index

			//find out if the nodes for this fiber have already been created
			it = find(id2vert.begin(), id2vert.end(), i[0]);	//look for the first node
			if (it == id2vert.end()) {							//if i[0] hasn't already been used
				vertex new_vertex = new_edge[0];				//create a new vertex, assign it a position
				new_vertex.e[0].push_back(E.size());			//add the current edge as outgoing
				new_edge.v[0] = V.size();						//add the new edge to the edge
				V.push_back(new_vertex);						//add the new vertex to the vertex list
				id2vert.push_back(i[0]);						//add the ID to the ID->vertex conversion list
			}
			else {									//if the vertex already exists
				it_idx = std::distance(id2vert.begin(), it);
				V[it_idx].e[0].push_back(E.size());				//add the current edge as outgoing
				new_edge.v[0] = it_idx;
			}

			it = find(id2vert.begin(), id2vert.end(), i[1]);	//look for the second ID
			if (it == id2vert.end()) {							//if i[1] hasn't already been used
				vertex new_vertex = new_edge[I - 1];			//create a new vertex, assign it a position
				new_vertex.e[1].push_back(E.size());			//add the current edge as incoming
				new_edge.v[1] = V.size();                       //add the new vertex to the edge
				V.push_back(new_vertex);						//add the new vertex to the vertex list
				id2vert.push_back(i[1]);						//add the ID to the ID->vertex conversion list
			}
			else {									//if the vertex already exists
				it_idx = std::distance(id2vert.begin(), it);
				V[it_idx].e[1].push_back(E.size());				//add the current edge as incoming
				new_edge.v[1] = it_idx;
			}

			E.push_back(new_edge);								//push the edge to the list
		}
	}

	/// Get adjacency matrix of the network
	std::vector< typename std::vector<int> > get_adj_mat() {
		
		unsigned n = V.size();		// get the number of vertices in the networks

		std::vector< typename std::vector<int> > result(n, std::vector<int>(n, 0));	// initialize every entry in the matrix to be 0
		result.resize(n);			// resize rows
		for (unsigned i = 0; i < n; i++)
			result[i].resize(n);	// resize columns
		
		for (unsigned i = 0; i < n; i++) {			// for every vertex
			unsigned num_out = V[i].e[0].size();	// number of outgoing edges of current vertex
			if (num_out != 0) {
				for (unsigned j = 0; j < num_out; j++) {
					int edge_idx = V[i].e[0][j];		// get the jth out-going edge index of current vertex
					int vertex_idx = E[edge_idx].v[1];	// get the ending vertex of specific out-going edge
					result[i][vertex_idx] = 1;			// can simply set to 1 if it is simple-graph
					result[vertex_idx][i] = 1;			// symmetric
				}
			}
		}

		return result;
	}

	/// Output the network as a string
	std::string str(){

		std::stringstream ss;
		ss<<"Nodes ("<<V.size()<<")--------"<<std::endl;
		for(unsigned int v = 0; v < V.size(); v++){
			ss<<"\t"<<v<<V[v].str()<<std::endl;
		}

		ss<<"Edges ("<<E.size()<<")--------"<<std::endl;
		for(unsigned e = 0; e < E.size(); e++){
			ss<<"\t"<<e<<E[e].str()<<std::endl;
		}

		return ss.str();
	}

	/// This function resamples all fibers in a network given a desired minimum spacing
	/// @param spacing is the minimum distance between two points on the network
	stim::network<T> resample(T spacing){
		stim::network<T> n;								//create a new network that will be an exact copy, with resampled fibers
		n.V = V;									//copy all vertices
		//n.NT = NT;										//copy all the neuronal type information
		n.E.resize(edges());								//allocate space for the edge list

		//copy all fibers, resampling them in the process
		for(unsigned e = 0; e < edges(); e++){						//for each edge in the edge list
			n.E[e] = E[e].resample(spacing);					//resample the edge and copy it to the new network
		}

		return n;							              	//return the resampled network
	}

	/// Calculate the total number of points on all edges.
	unsigned total_points(){
		unsigned n = 0;
		for(unsigned e = 0; e < E.size(); e++)
			n += E[e].size();
		return n;
	}

	//Copy the point cloud representing the centerline for the network into an array
	void centerline_cloud(T* dst) {
		size_t p;										//stores the current edge point
		size_t P;										//stores the number of points in an edge
		size_t i = 0;									//index into the output array of points
		for (size_t e = 0; e < E.size(); e++) {			//for each edge in the network
			P = E[e].size();							//get the number of points in this edge
			for (p = 0; p < P; p++) {
				dst[i * 3 + 0] = E[e][p][0];		
				dst[i * 3 + 1] = E[e][p][1];
				dst[i * 3 + 2] = E[e][p][2];
				i++;
			}
		}
	}

    // convert vec3 to array
	void stim2array(float *a, stim::vec3<T> b){
		a[0] = b[0];
		a[1] = b[1];
		a[2] = b[2];
	}

	// convert vec3 to array in bunch
	void edge2array(T* a, edge b){
		size_t n = b.size();
		for(size_t i = 0; i < n; i++){
			a[i * 3 + 0] = b[i][0];
			a[i * 3 + 1] = b[i][1];
			a[i * 3 + 2] = b[i][2];	 
		}
	}

	// get list of metric
	std::vector<T> metric() {
		std::vector<T> result;
		T m;
		for (size_t e = 0; e < E.size(); e++) {
			for (size_t p = 0; p < E[e].size(); p++) {
				m = E[e].r(p);
				result.push_back(m);
			}
		}
		return result;
	}

	/// Calculate the average magnitude across the entire network.
	/// @param m is the magnitude value to use. The default is 0 (usually radius).
	T average(unsigned m = 0){

		T M, L;										//allocate space for the total magnitude and length
		M = L = 0;									//initialize both the initial magnitude and length to zero
		for(unsigned e = 0; e < E.size(); e++){						//for each edge in the network
			M += E[e].integrate();							//get the integrated magnitude
			L += E[e].length();							//get the edge length
		}

		return M / L;									//return the average magnitude
	}

	/// This function compares two networks and returns the percentage of the current network that is missing from A.
	/// @param A is the network to compare to - the field is generated for A
	/// @param sigma is the user-defined tolerance value - smaller values provide a stricter comparison
	stim::network<T> compare(stim::network<T> A, float sigma, int device = -1){

		stim::network<T> R;										//generate a network storing the result of the comparison
		R = (*this);											//initialize the result with the current network

		T *c;						                 			// centerline (array of double pointers) - points on kdtree must be double
		size_t n_data = A.total_points();          				// set the number of points
		c = (T*) malloc(sizeof(T) * n_data * 3);				// allocate an array to store all points in the data set				

		unsigned t = 0;
		for(unsigned e = 0; e < A.E.size(); e++){				//for each edge in the network
			for(unsigned p = 0; p < A.E[e].size(); p++){		//for each point in the edge
				for(unsigned d = 0; d < 3; d++){				//for each coordinate

					c[t * 3 + d] = A.E[e][p][d];				//copy the point into the array c
				}
				t++;
			}
		}

		//generate a KD-tree for network A
		size_t MaxTreeLevels = 3;								// max tree level
		
#ifdef __CUDACC__
		cudaSetDevice(device);
		int current_device;
		if (cudaGetDevice(&current_device) == device) {
			std::cout << "Using CUDA device " << device << " for calculations..." << std::endl;
		}
		stim::kdtree<T, 3> kdt;								// initialize a pointer to a kd tree

		kdt.create(c, n_data, MaxTreeLevels);				// build a KD tree

		for(unsigned e = 0; e < R.E.size(); e++){					//for each edge in A
			//R.E[e].add_mag(0);							//add a new magnitude for the metric
			//size_t errormag_id = R.E[e].nmags() - 1;		//get the id for the new magnitude
			
			size_t n = R.E[e].size();						// the number of points in current edge
			T* queryPt = new T[3 * n];
			T* m1 = new T[n];
			T* dists = new T[n];
			size_t* nnIdx = new size_t[n];

			T* d_dists;										
			T* d_m1;										
			cudaMalloc((void**)&d_dists, n * sizeof(T));
			cudaMalloc((void**)&d_m1, n * sizeof(T));

			edge2array(queryPt, R.E[e]);
			kdt.search(queryPt, n, nnIdx, dists);		

			cudaMemcpy(d_dists, dists, n * sizeof(T), cudaMemcpyHostToDevice);					// copy dists from host to device

			// configuration parameters
			size_t threads = (1024>n)?n:1024;
			size_t blocks = n/threads + (n%threads)?1:0;

			find_metric_parallel<<<blocks, threads>>>(d_m1, n, d_dists, sigma);					//calculate the metric value based on the distance

			cudaMemcpy(m1, d_m1, n * sizeof(T), cudaMemcpyDeviceToHost);

			for(unsigned p = 0; p < n; p++){
				R.E[e].set_r(p, m1[p]);
			}

			//d_set_mag<<<blocks, threads>>>(R.E[e].M, errormag_id, n, m1);
		}

#else
		stim::kdtree<T, 3> kdt;
		kdt.create(c, n_data, MaxTreeLevels);
	
		for(unsigned e = 0; e < R.E.size(); e++){			//for each edge in A

			size_t n = R.E[e].size();						// the number of points in current edge
			T* query = new T[3 * n];
			T* m1 = new T[n];
			T* dists = new T[n];
			size_t* nnIdx = new size_t[n];

			edge2array(query, R.E[e]);

			kdt.cpu_search(query, n, nnIdx, dists);			//find the distance between A and the current network

			for(unsigned p = 0; p < R.E[e].size(); p++){
				m1[p] = 1.0f - gaussianFunction((T)dists[p], sigma);	//calculate the metric value based on the distance
				R.E[e].set_r(p, m1[p]);					//set the error for the second point in the segment
			}
		}
#endif
		return R;		//return the resulting network
	}

	/// This function compares two networks and split the current one according to the nearest neighbor of each point in each edge
	/// @param A is the network to split
	/// @param B is the corresponding mapping network
	/// @param sigma is the user-defined tolerance value - smaller values provide a stricter comparison
	/// @param device is the device that user want to use
	void split(stim::network<T> A, stim::network<T> B, float sigma, int device, float threshold){

		T *c;						                 	
		size_t n_data = B.total_points();          				
		c = (T*) malloc(sizeof(T) * n_data * 3); 				

		size_t NB = B.E.size();								// the number of edges in B
		unsigned t = 0;
		for(unsigned e = 0; e < NB; e++){					// for every edge in B			
			for(unsigned p = 0; p < B.E[e].size(); p++){	// for every points in B.E[e]
				for(unsigned d = 0; d < 3; d++){				

					c[t * 3 + d] = B.E[e][p][d];			// convert to array
				}
				t++;
			}
		}
		size_t MaxTreeLevels = 3;							// max tree level

#ifdef __CUDACC__
		cudaSetDevice(device);
		stim::kdtree<T, 3> kdt;								// initialize a pointer to a kd tree
	
		//compare each point in the current network to the field produced by A
		kdt.create(c, n_data, MaxTreeLevels);				// build a KD tree

		std::vector<std::vector<unsigned> > relation;		// the relationship between GT and T corresponding to NN
		relation.resize(A.E.size());										

		for(unsigned e = 0; e < A.E.size(); e++){			//for each edge in A
			//A.E[e].add_mag(0);								//add a new magnitude for the metric
			//size_t errormag_id = A.E[e].nmags() - 1;
			
			size_t n = A.E[e].size();						// the number of edges in A

			T* queryPt = new T[3 * n];							// set of all the points in current edge
			T* m1 = new T[n];								// array of metrics for every point in current edge
			T* dists = new T[n];							// store the distances for every point in current edge
			size_t* nnIdx = new size_t[n];					// store the indices for every point in current edge
			
			// define pointers in device
			T* d_dists;														
			T* d_m1;
			size_t* d_nnIdx;

			// allocate memory for defined pointers
			cudaMalloc((void**)&d_dists, n * sizeof(T));
			cudaMalloc((void**)&d_m1, n * sizeof(T));
			cudaMalloc((void**)&d_nnIdx, n * sizeof(size_t));

			edge2array(queryPt, A.E[e]);						// convert edge to array
			kdt.search(queryPt, n, nnIdx, dists);				// search the tree to find the NN for every point in current edge

			cudaMemcpy(d_dists, dists, n * sizeof(T), cudaMemcpyHostToDevice);					// copy dists from host to device
			cudaMemcpy(d_nnIdx, nnIdx, n * sizeof(size_t), cudaMemcpyHostToDevice);				// copy Idx from host to device

			// configuration parameters
			size_t threads = (1024>n)?n:1024;													// test to see whether the number of point in current edge is more than 1024
			size_t blocks = n/threads + (n%threads)?1:0;

			find_metric_parallel<<<blocks, threads>>>(d_m1, n, d_dists, sigma);								// calculate the metrics in parallel

			cudaMemcpy(m1, d_m1, n * sizeof(T), cudaMemcpyDeviceToHost);

			for(unsigned p = 0; p < n; p++){
				A.E[e].set_r(p, m1[p]);											// set the error(radius) value to every point in current edge
			}

			relation[e].resize(n);																// resize every edge relation size

			unsigned* d_relation;
			cudaMalloc((void**)&d_relation, n * sizeof(unsigned));								// allocate memory

			std::vector<size_t> edge_point_num(NB);												// %th element is the number of points that %th edge has
			for(unsigned ee = 0; ee < NB; ee++)
				edge_point_num[ee] = B.E[ee].size();

			size_t* d_edge_point_num;
			cudaMalloc((void**)&d_edge_point_num, NB * sizeof(size_t));
			cudaMemcpy(d_edge_point_num, &edge_point_num[0], NB * sizeof(size_t), cudaMemcpyHostToDevice);

			find_edge_index_parallel<<<blocks, threads>>>(d_nnIdx, n, d_relation, d_edge_point_num, NB);			// find the edge corresponding to the array index in parallel

			cudaMemcpy(&relation[e][0], d_relation, n * sizeof(unsigned), cudaMemcpyDeviceToHost);	//copy relationship from device to host
		}
#else
		stim::kdtree<T, 3> kdt;
		kdt.create(c, n_data, MaxTreeLevels);
	
		std::vector<std::vector<unsigned>> relation;		// the mapping relationship between two networks
		relation.resize(A.E.size());										
		for(unsigned i = 0; i < A.E.size(); i++)
			relation[i].resize(A.E[i].size());

		std::vector<size_t> edge_point_num(NB);				//%th element is the number of points that %th edge has
		for(unsigned ee = 0; ee < NB; ee++)
			edge_point_num[ee] = B.E[ee].size();

		for(unsigned e = 0; e < A.E.size(); e++){			//for each edge in A
			
			size_t n = A.E[e].size();						//the number of edges in A

			T* queryPt = new T[3 * n];
			T* m1 = new T[n];
			T* dists = new T[n];							//store the distances
			size_t* nnIdx = new size_t[n];					//store the indices
			
			edge2array(queryPt, A.E[e]);
			kdt.search(queryPt, n, nnIdx, dists);		

			for(unsigned p = 0; p < A.E[e].size(); p++){
				m1[p] = 1.0f - gaussianFunction((T)dists[p], sigma);	//calculate the metric value based on the distance
				A.E[e].set_r(p, m1[p]);									//set the error for the second point in the segment
				
				unsigned id = 0;																	//mapping edge's idx
				size_t num = 0;																		//total number of points before #th edge
				for(unsigned i = 0; i < NB; i++){
					num += B.E[i].size();
					if(nnIdx[p] < num){																//find the edge it belongs to
						relation[e][p] = id;
						break;
					}
					id++;																			//current edge won't be the one, move to next edge
				}
			}
		}
#endif
		E = A.E;
		V = A.V;

		unsigned int id = 0;									// split value								
		for(unsigned e = 0; e < E.size(); e++){					// for every edge
			for(unsigned p = 0; p < E[e].size() - 1; p++){		// for every point in each edge
				int t = (int)(E[e].length() / sigma * 2);
				if (t <= 20)
					threshold_fac = E[e].size();
				else
					threshold_fac = (E[e].length() / sigma * 2)/10;
				if(relation[e][p] != relation[e][p + 1]){		// find the nearest edge changing point
					id = p + 1;									// if there is no change in NN
					if(id < threshold_fac || (E[e].size() - id) < threshold_fac)				
						id = E[e].size() - 1;																			// extreme situation is not acceptable
					else
						break;
				}
				if(p == E[e].size() - 2)																// if there is no splitting index, set the id to the last point index of current edge
					id = E[e].size() - 1;
			}
			//unsigned errormag_id = E[e].nmags() - 1;
			T G = 0;																					// test to see whether it has its nearest neighbor
			for(unsigned i = 0; i < E[e].size(); i++)
				G += E[e].r(i);																			// won't split special edges
			if(G / E[e].size() > threshold)															// should based on the color map
				id = E[e].size() - 1;																	// set split idx to outgoing direction vertex

			std::vector<edge> tmpe;
			tmpe.resize(2);
			tmpe = E[e].split(id);
			vertex tmpv = stim::vec3<T>(-1, -1, 0);														// store the split point as vertex
			if(tmpe.size() == 2){
				relation.resize(relation.size() + 1);
				for(unsigned d = id; d < E[e].size(); d++)
					relation[relation.size() - 1].push_back(relation[e][d]);
				tmpe[0].v[0] = E[e].v[0];																// begining vertex of first half edge -> original begining vertex
				tmpe[1].v[1] = E[e].v[1];																// ending vertex of second half edge -> original ending vertex
				tmpv = E[e][id];
				V.push_back(tmpv);
				tmpe[0].v[1] = (unsigned)V.size() - 1;													// ending vertex of first half edge -> new vertex
				tmpe[1].v[0] = (unsigned)V.size() - 1;													// begining vertex of second half edge -> new vertex
				edge tmp(E[e]);
				E[e] = tmpe[0];																			// replace original edge by first half edge
				E.push_back(tmpe[1]);																	// push second half edge to the last
				V[V.size() - 1].e[1].push_back(e);														// push first half edge to the incoming of new vertex
				V[V.size() - 1].e[0].push_back((unsigned)E.size() - 1);									// push second half edge to the outgoing of new vertex
				for(unsigned i = 0; i < V[tmp.v[1]].e[1].size(); i++)									// find the incoming edge of original ending vertex
					if(V[tmp.v[1]].e[1][i] == e)
						V[tmp.v[1]].e[1][i] = (unsigned)E.size() - 1;									// set to new edge
			}
		}
	}

	/// This function compares two splitted networks and yields a mapping relationship between them according to NN
	/// @param B is the network that the current network is going to map to
	/// @param C is the mapping relationship: C[e1] = _e1 means e1 edge in current network is mapping the _e1 edge in B
	/// @param device is the device that user want to use
	void mapping(stim::network<T> B, std::vector<unsigned> &C, int device, float threshold){
		stim::network<T> A;								//generate a network storing the result of the comparison
		A = (*this);

		size_t n = A.E.size();							// the number of edges in A
		size_t NB = B.E.size();							// the number of edges in B

		C.resize(A.E.size());	

		T *c;						                 	// centerline (array of double pointers) - points on kdtree must be double
		size_t n_data = B.total_points();          		// set the number of points
		c = (T*) malloc(sizeof(T) * n_data * 3); 				

		unsigned t = 0;
		for(unsigned e = 0; e < NB; e++){					// for each edge in the network
			for(unsigned p = 0; p < B.E[e].size(); p++){	// for each point in the edge
				for(unsigned d = 0; d < 3; d++){			// for each coordinate

					c[t * 3 + d] = B.E[e][p][d];
				}
				t++;
			}
		}

		//generate a KD-tree for network A
		//float metric = 0.0;                               		// initialize metric to be returned after comparing the network
		size_t MaxTreeLevels = 3;									// max tree level
		
#ifdef __CUDACC__
		cudaSetDevice(device);
		stim::kdtree<T, 3> kdt;								// initialize a pointer to a kd tree
	
		kdt.create(c, n_data, MaxTreeLevels);				// build a KD tree

		for(unsigned e = 0; e < n; e++){					//for each edge in A
			//size_t errormag_id = A.E[e].nmags() - 1;		//get the id for the new magnitude
			
			//pre-judge to get rid of impossibly mapping edges
			T M = 0;
			for(unsigned p = 0; p < A.E[e].size(); p++)
				M += A.E[e].r(p);
			M = M / A.E[e].size();
			if(M > threshold)
				C[e] = (unsigned)-1;						//set the nearest edge of impossibly mapping edges to maximum of unsigned
			else{
				T* queryPt = new T[3];
				T* dists = new T[1];
				size_t* nnIdx = new size_t[1];

				stim2array(queryPt, A.E[e][A.E[e].size()/2]);
				kdt.search(queryPt, 1, nnIdx, dists);
				
				unsigned id = 0;							//mapping edge's idx
				size_t num = 0;								//total number of points before #th edge
				for(unsigned i = 0; i < NB; i++){
					num += B.E[i].size();
					if(nnIdx[0] < num){
						C[e] = id;
						break;
					}
					id++;
				}
			}
		}
#else
		stim::kdtree<T, 3> kdt;
		kdt.create(c, n_data, MaxTreeLevels);
		T *dists = new T[1];								// near neighbor distances
		size_t *nnIdx = new size_t[1];						// near neighbor indices // allocate near neigh indices

		stim::vec3<T> p0, p1;
		T* queryPt = new T[3];

		for(unsigned int e = 0; e < R.E.size(); e++){			// for each edge in A
			T M;											// the sum of metrics of current edge
			for(unsigned p = 0; p < R.E[e].size(); p++)
				M += A.E[e].r(p);
			M = M / A.E[e].size();
			if(M > threshold)								
				C[e] = (unsigned)-1;
			else{											// if it should have corresponding edge in B, then...
				p1 = R.E[e][R.E[e].size()/2];							
				stim2array(queryPt, p1);
				kdt.cpu_search(queryPt, 1, nnIdx, dists);	// search the tree		
				
				unsigned id = 0;							//mapping edge's idx
				size_t num = 0;								//total number of points before #th edge
				for(unsigned i = 0; i < NB; i++){
					num += B.E[i].size();
					if(nnIdx[0] < num){
						C[e] = id;
						break;
					}
					id++;
				}
			}
		}
#endif
	}

	/// Returns the number of magnitude values stored in each edge. This should be uniform across the network.
	//unsigned nmags(){
	//	return E[0].nmags();
	//}
	// split a string in text by the character sep
	stim::vec<T> split(std::string &text, char sep) 
	{
		stim::vec<T> tokens;
		std::size_t start = 0, end = 0;
		while ((end = text.find(sep, start)) != std::string::npos) {
		tokens.push_back(atof(text.substr(start, end - start).c_str()));
		start = end + 1;
		}
		tokens.push_back(atof(text.substr(start).c_str()));
		return tokens;
	}
	// load a network in text file to a network class
	void load_txt(std::string filename)
	{
		std::vector <std::string> file_contents;
		std::ifstream file(filename.c_str());
		std::string line;
		std::vector<unsigned> id2vert;	//this list stores the vertex ID associated with each network vertex
		//for each line in the text file, store them as strings in file_contents
		while (std::getline(file, line))
		{
			std::stringstream ss(line);
			file_contents.push_back(ss.str());
		}
		unsigned int numEdges = atoi(file_contents[0].c_str()); //number of edges in the network
		unsigned int I = atoi(file_contents[1].c_str()) ;		//calculate the number of points3d on the first edge
		unsigned int count = 1; unsigned int k = 2; // count is global counter through the file contents, k is for the vertices on the edges
		// for each edge in the network.
		for (unsigned int i = 0; i < numEdges; i ++ )
		{
			// pre allocate a position vector p with number of points3d on the edge p
			std::vector< stim::vec<T> > p(0, I);
			// for each point on the nth edge
		  for (unsigned int j = k; j < I + k; j++)
		  {
			 // split the points3d of floats with separator space and form a float3 position vector out of them
			  p.push_back(split(file_contents[j], ' '));
		  }
			count += p.size() + 1; // increment count to point at the next edge in the network
			I = atoi(file_contents[count].c_str()); // read in the points3d at the next edge and convert it to an integer
			k = count + 1;
			edge new_edge = p; // create an edge with a vector of points3d  on the edge
			E.push_back(new_edge); // push the edge into the network
		}
		unsigned int numVertices = atoi(file_contents[count].c_str()); // this line in the text file gives the number of distinct vertices
		count = count + 1; // this line of text file gives the first verrtex
		// push each vertex into V
		for (unsigned int i = 0; i < numVertices; i ++)
		{
			vertex new_vertex = split(file_contents[count], ' ');
			V.push_back(new_vertex);
			count += atoi(file_contents[count + 1].c_str()) + 2; // Skip number of edge ids + 2 to point to the next vertex
		}
	} // end load_txt function

	// strTxt returns a string of edges
	std::string
	strTxt(std::vector< stim::vec<T> > p)
	{
		std::stringstream ss;
		std::stringstream oss;
		for(unsigned int i = 0; i < p.size(); i++){
			ss.str(std::string());
			for(unsigned int d = 0; d < 3; d++){
				ss<<p[i][d];
			}
			ss << "\n";
		}
		return ss.str();
	}
	// removes specified character from string
	void removeCharsFromString(std::string &str, char* charsToRemove ) {
	   for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
		  str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
	   }
	}
	//exports network to txt file
	void
	to_txt(std::string filename)
	{
		std::ofstream ofs(filename.c_str(), std::ofstream::out | std::ofstream::app);
		//int num;
		ofs << (E.size()).str() << "\n";
		for(unsigned int i = 0; i < E.size(); i++)
		{
			 std::string str;
			 ofs << (E[i].size()).str() << "\n";
			 str = E[i].strTxt();
             ofs << str << "\n";
		} 
		for(int i = 0; i < V.size(); i++)
		{
			std::string str;
			str = V[i].str();
			char temp[4] = "[],";
			removeCharsFromString(str, temp);
			ofs << str << "\n";
		}
		ofs.close();
	}
};		//end stim::network class
};		//end stim namespace
#endif
