#ifndef STIM_CENTERLINE_H
#define STIM_CENTERLINE_H

#include <vector>
#include <stim/math/vec3.h>
#include <stim/structures/kdtree.cuh>

namespace stim{

/**	This class stores information about a single fiber represented as a set of geometric points
 *	between two branch or end points. This class is used as a fundamental component of the stim::network
 *	class to describe an interconnected (often biological) network.
 */
template<typename T>
class centerline : public std::vector< stim::vec3<T> >{

protected:

	std::vector<T> L;										//stores the integrated length along the fiber (used for parameterization)

	///Return the normalized direction vector at point i (average of the incoming and outgoing directions)
	vec3<T> d(size_t i) {
		if (size() <= 1) return vec3<T>(0, 0, 0);						//if there is insufficient information to calculate the direction, return a null vector
		if (size() == 2) return (at(1) - at(0)).norm();					//if there are only two points, the direction vector at both is the direction of the line segment
		if (i == 0) return (at(1) - at(0)).norm();						//the first direction vector is oriented towards the first line segment
		if (i == size() - 1) return (at(size() - 1) - at(size() - 2)).norm();	//the last direction vector is oriented towards the last line segment

		//all other direction vectors are the average direction of the two joined line segments
		vec3<T> a = at(i) - at(i - 1);
		vec3<T> b = at(i + 1) - at(i);
		vec3<T> ab = a.norm() + b.norm();
		return ab.norm();
	}

	//initializes the integrated length vector to make parameterization easier, starting with index idx (all previous indices are assumed to be correct)
	void update_L(size_t start = 0) {
		L.resize(size());									//allocate space for the L array
		if (start == 0) {
			L[0] = 0;											//initialize the length value for the first point to zero (0)
			start++;
		}

		stim::vec3<T> d;
		for (size_t i = start; i < size(); i++) {		//for each line segment in the centerline
			d = at(i) - at(i - 1);
			L[i] = L[i - 1] + d.len();				//calculate the running length total
		}
	}

	void init() {
		if (size() == 0) return;								//return if there aren't any points
		update_L();
	}

	/// Returns a stim::vec representing the point at index i

	/// @param i is an index of the desired centerline point
	stim::vec<T> get_vec(unsigned i){
		return std::vector< stim::vec3<T> >::at(i);
	}

	///finds the index of the point closest to the length l on the lower bound.
	///binary search.
	size_t findIdx(T l) {
		for (size_t i = 1; i < L.size(); i++) {				//for each point in the centerline
			if (L[i] > l) return i - 1;						//if we have passed the desired length value, return i-1
		}
		return L.size() - 1;
		/*size_t i = L.size() / 2;
		size_t max = L.size() - 1;
		size_t min = 0;
		while (i < L.size() - 1){
			if (l < L[i]) {
				max = i;
				i = min + (max - min) / 2;
			}
			else if (L[i] <= l && L[i + 1] >= l) {
				break;
			}
			else {
				min = i;
				i = min + (max - min) / 2;
			}
		}
		return i;*/
	}

	///Returns a position vector at the given length into the fiber (based on the pvalue).
	///Interpolates the radius along the line.
	///@param l: the location of the in the cylinder.
	///@param idx: integer location of the point closest to l but prior to it.
	stim::vec3<T> p(T l, int idx) {
		T rat = (l - L[idx]) / (L[idx + 1] - L[idx]);
		stim::vec3<T> v1 = at(idx);
		stim::vec3<T> v2 = at(idx + 1);
		return(v1 + (v2 - v1)*rat);
	}


public:

	using std::vector< stim::vec3<T> >::at;
	using std::vector< stim::vec3<T> >::size;

	centerline() : std::vector< stim::vec3<T> >() {
		init();
	}
	centerline(size_t n) : std::vector< stim::vec3<T> >(n){
		init();
	}
	centerline(std::vector<stim::vec3<T> > pos) :
		std::vector<stim::vec3<T> > (pos)
	{
		init();
	}
	
	//overload the push_back function to update the length vector
	void push_back(stim::vec3<T> p) {
		std::vector< stim::vec3<T> >::push_back(p);
		update_L(size() - 1);
	}

	///Returns a position vector at the given p-value (p value ranges from 0 to 1).
	///interpolates the position along the line.
	///@param pvalue: the location of the in the cylinder, from 0 (beginning to 1).
	stim::vec3<T> p(T pvalue) {
		if (pvalue <= 0.0) return at(0);			//return the first element
		if (pvalue >= 1.0) return back();			//return the last element

		T l = pvalue*L[L.size() - 1];
		int idx = findIdx(l);
		return p(l, idx);
	}

	///Update centerline internal parameters (currently the L vector)
	void update() {
		init();
	}
	///Return the length of the entire centerline
	T length() {
		return L.back();
	}


	/// stitch two centerlines
	///@param c1, c2: two centerlines
	///@param sigma: sample rate
	static std::vector< stim::centerline<T> > stitch(stim::centerline<T> c1, stim::centerline<T> c2 = stim::centerline<T>()) {
		
		std::vector< stim::centerline<T> > result;
		stim::centerline<T> new_centerline;
		stim::vec3<T> new_vertex;

		// if only one centerline, stitch itself!
		if (c2.size() == 0) {
			size_t num = c1.size();
			size_t id = 100000;							// store the downsample start position
			T threshold;
			if (num < 4) {								// if the number of vertex is less than 4, do nothing
				result.push_back(c1);
				return result;
			}
			else {
				// test geometry start vertex
				stim::vec3<T> v1 = c1[1] - c1[0];		// vector from c1[0] to c1[1]
				for (size_t p = 2; p < num; p++) {		// 90° standard???
					stim::vec3<T> v2 = c1[p] - c1[0];
					float cosine = v2.dot(v1);
					if (cosine < 0) {
						id = p;
						threshold = v2.len();
						break;
					}
				}
				if (id != 100000) {						// find a downsample position on the centerline
					T* c;
					c = (T*)malloc(sizeof(T) * (num - id) * 3);
					for (size_t p = id; p < num; p++) {
						for (size_t d = 0; d < 3; d++) {
							c[p * 3 + d] = c1[p][d];
						}
					}
					stim::kdtree<T, 3> kdt;
					kdt.create(c, num - id, 5);			// create tree

					T* query = (T*)malloc(sizeof(T) * 3);
					for (size_t d = 0; d < 3; d++)
						query[d] = c1[0][d];
					size_t index;
					T dist;

					kdt.search(query, 1, &index, &dist);

					free(query);
					free(c);

					if (dist > threshold) {
						result.push_back(c1);
					}
					else {
						// the loop part
						new_vertex = c1[index];
						new_centerline.push_back(new_vertex);
						for (size_t p = 0; p < index + 1; p++) {
							new_vertex = c1[p];
							new_centerline.push_back(new_vertex);
						}
						result.push_back(new_centerline);
						new_centerline.clear();

						// the tail part
						for (size_t p = index; p < num; p++) {
							new_vertex = c1[p];
							new_centerline.push_back(new_vertex);
						}
						result.push_back(new_centerline);
					}
				}
				else {	// there is one potential problem that two positions have to be stitched
						// test geometry end vertex
					stim::vec3<T> v1 = c1[num - 2] - c1[num - 1];
					for (size_t p = num - 2; p > 0; p--) {		// 90° standard
						stim::vec3<T> v2 = c1[p - 1] - c1[num - 1];
						float cosine = v2.dot(v1);
						if (cosine < 0) {
							id = p;
							threshold = v2.len();
							break;
						}
					}
					if (id != 100000) {						// find a downsample position
						T* c;
						c = (T*)malloc(sizeof(T) * (id + 1) * 3);
						for (size_t p = 0; p < id + 1; p++) {
							for (size_t d = 0; d < 3; d++) {
								c[p * 3 + d] = c1[p][d];
							}
						}
						stim::kdtree<T, 3> kdt;
						kdt.create(c, id + 1, 5);				// create tree

						T* query = (T*)malloc(sizeof(T) * 1 * 3);
						for (size_t d = 0; d < 3; d++)
							query[d] = c1[num - 1][d];
						size_t index;
						T dist;

						kdt.search(query, 1, &index, &dist);

						free(query);
						free(c);

						if (dist > threshold) {
							result.push_back(c1);
						}
						else {
							// the tail part
							for (size_t p = 0; p < index + 1; p++) {
								new_vertex = c1[p];
								new_centerline.push_back(new_vertex);
							}
							result.push_back(new_centerline);
							new_centerline.clear();

							// the loop part
							for (size_t p = index; p < num; p++) {
								new_vertex = c1[p];
								new_centerline.push_back(new_vertex);
							}
							new_vertex = c1[index];
							new_centerline.push_back(new_vertex);
							result.push_back(new_centerline);
						}
					}
					else {	// no stitch position
						result.push_back(c1);
					}
				}
			}
		}


		// two centerlines
		else {
			// find stitch position based on nearest neighbors												
			size_t num1 = c1.size();
			T* c = (T*)malloc(sizeof(T) * num1 * 3);		// c1 as reference point
			for (size_t p = 0; p < num1; p++)				// centerline to array
				for (size_t d = 0; d < 3; d++)				// because right now my kdtree code is a relatively close code, it has its own structure
					c[p * 3 + d] = c1[p][d];				// I will merge it into stimlib totally in the near future

			stim::kdtree<T, 3> kdt;							// kdtree object
			kdt.create(c, num1, 5);							// create tree

			size_t num2 = c2.size();						
			T* query = (T*)malloc(sizeof(T) * num2 * 3);	// c2 as query point
			for (size_t p = 0; p < num2; p++) {
				for (size_t d = 0; d < 3; d++) {
					query[p * 3 + d] = c2[p][d];
				}
			}
			std::vector<size_t> index(num2);
			std::vector<T> dist(num2);

			kdt.search(query, num2, &index[0], &dist[0]);	// find the nearest neighbors in c1 for c2

			// clear up
			free(query);
			free(c);

			// find the average vertex distance of one centerline
			T sigma1 = 0;
			T sigma2 = 0;
			for (size_t p = 0; p < num1 - 1; p++) 
				sigma1 += (c1[p] - c1[p + 1]).len();
			for (size_t p = 0; p < num2 - 1; p++)
				sigma2 += (c2[p] - c2[p + 1]).len();
			sigma1 /= (num1 - 1);
			sigma2 /= (num2 - 1);
			float threshold = 4 * (sigma1 + sigma2) / 2;			// better way to do this?

			T min_d = *std::min_element(dist.begin(), dist.end());	// find the minimum distance between c1 and c2

			if (min_d > threshold) {								// if the minimum distance is too large
				result.push_back(c1);
				result.push_back(c2);

#ifdef DEBUG
				std::cout << "The distance between these two centerlines is too large" << std::endl;
#endif
			}
			else {
			//	auto smallest = std::min_element(dist.begin(), dist.end());
				unsigned int smallest = std::min_element(dist.begin(), dist.end());
			//	auto i = std::distance(dist.begin(), smallest);		// find the index of min-distance in distance list
				unsigned int i = std::distance(dist.begin(), smallest);		// find the index of min-distance in distance list

				// stitch position in c1 and c2
				int id1 = index[i];
				int id2 = i;

				// actually there are two cases
				// first one inacceptable
				// second one acceptable
				if (id1 != 0 && id1 != num1 - 1 && id2 != 0 && id2 != num2 - 1) {		// only stitch one end vertex to another centerline
					result.push_back(c1);
					result.push_back(c2);
				}
				else {
					if (id1 == 0 || id1 == num1 - 1) {			// if the stitch vertex is the first or last vertex of c1
						// for c2, consider two cases(one degenerate case)
						if (id2 == 0 || id2 == num2 - 1) {		// case 1, if stitch position is also on the end of c2
							// we have to decide which centerline get a new vertex, based on direction
							// for c1, computer the direction change angle
							stim::vec3<T> v1, v2;
							float alpha1, alpha2;				// direction change angle
							if (id1 == 0)
								v1 = (c1[1] - c1[0]).norm();
							else
								v1 = (c1[num1 - 2] - c1[num1 - 1]).norm();
							v2 = (c2[id2] - c1[id1]).norm();
							alpha1 = v1.dot(v2);
							if (id2 == 0)
								v1 = (c2[1] - c2[0]).norm();
							else
								v1 = (c2[num2 - 2] - c2[num2 - 1]).norm();
							v2 = (c1[id1] - c2[id2]).norm();
							alpha2 = v1.dot(v2);
							if (abs(alpha1) > abs(alpha2)) {					// add the vertex to c1 in order to get smooth connection
								// push back c1
								if (id1 == 0) {									// keep geometry information
									new_vertex = c2[id2];
									new_centerline.push_back(new_vertex);
									for (size_t p = 0; p < num1; p++) {			// stitch vertex on c2 -> geometry start vertex on c1 -> geometry end vertex on c1
										new_vertex = c1[p];
										new_centerline.push_back(new_vertex);
									}
								}
								else {
									for (size_t p = 0; p < num1; p++) {			// stitch vertex on c2 -> geometry end vertex on c1 -> geometry start vertex on c1
										new_vertex = c1[p];
										new_centerline.push_back(new_vertex);
									}
									new_vertex = c2[id2];
									new_centerline.push_back(new_vertex);
								}
								result.push_back(new_centerline);
								new_centerline.clear();

								// push back c2
								for (size_t p = 0; p < num2; p++) {
									new_vertex = c2[p];
									new_centerline.push_back(new_vertex);
								}
								result.push_back(new_centerline);
							}
							else {												// add the vertex to c2 in order to get smooth connection
								// push back c1
								for (size_t p = 0; p < num1; p++) {
									new_vertex = c1[p];
									new_centerline.push_back(new_vertex);
								}
								result.push_back(new_centerline);
								new_centerline.clear();

								// push back c2
								if (id2 == 0) {									// keep geometry information
									new_vertex = c1[id1];
									new_centerline.push_back(new_vertex);
									for (size_t p = 0; p < num2; p++) {			// stitch vertex on c2 -> geometry start vertex on c1 -> geometry end vertex on c1
										new_vertex = c2[p];
										new_centerline.push_back(new_vertex);
									}
								}
								else {
									for (size_t p = 0; p < num2; p++) {			// stitch vertex on c2 -> geometry end vertex on c1 -> geometry start vertex on c1
										new_vertex = c2[p];
										new_centerline.push_back(new_vertex);
									}
									new_vertex = c1[id1];
									new_centerline.push_back(new_vertex);
								}
								result.push_back(new_centerline);
							}
						}
						else {												// case 2, the stitch position is on c2
							// push back c1
							if (id1 == 0) {									// keep geometry information
								new_vertex = c2[id2];
								new_centerline.push_back(new_vertex);
								for (size_t p = 0; p < num1; p++) {			// stitch vertex on c2 -> geometry start vertex on c1 -> geometry end vertex on c1
									new_vertex = c1[p];
									new_centerline.push_back(new_vertex);
								}
							}
							else {
								for (size_t p = 0; p < num1; p++) {			// geometry end vertex on c1 -> geometry start vertex on c1 -> stitch vertex on c2
									new_vertex = c1[p];
									new_centerline.push_back(new_vertex);
								}
								new_vertex = c2[id2];
								new_centerline.push_back(new_vertex);
							}
							result.push_back(new_centerline);
							new_centerline.clear();

							// push back c2
							for (size_t p = 0; p < id2 + 1; p++) {			// first part
								new_vertex = c2[p];
								new_centerline.push_back(new_vertex);
							}
							result.push_back(new_centerline);
							new_centerline.clear();

							for (size_t p = id2; p < num2; p++) {			// second part
								new_vertex = c2[p];
								new_centerline.push_back(new_vertex);
							}
							result.push_back(new_centerline);
						}
					}
					else {							// if the stitch vertex is the first or last vertex of c2
						// push back c2
						if (id2 == 0) {										// keep geometry information
							new_vertex = c1[id1];
							new_centerline.push_back(new_vertex);
							for (size_t p = 0; p < num2; p++) {				// stitch vertex on c1 -> geometry start vertex on c2 -> geometry end vertex on c2
								new_vertex = c2[p];
								new_centerline.push_back(new_vertex);
							}
						}
						else {
							for (size_t p = 0; p < num2; p++) {				// geometry end vertex on c2 -> geometry start vertex on c2 -> stitch vertex on c1
								new_vertex = c2[p];
								new_centerline.push_back(new_vertex);
							}
							new_vertex = c1[id1];
							new_centerline.push_back(new_vertex);
							result.push_back(new_centerline);
							new_centerline.clear();

							// push back c1
							for (size_t p = 0; p < id1 + 1; p++) {			// first part
								new_vertex = c1[p];
								new_centerline.push_back(new_vertex);
							}
							result.push_back(new_centerline);
							new_centerline.clear();

							for (size_t p = id1; p < num1; p++) {			// second part
								new_vertex = c1[p];
								new_centerline.push_back(new_vertex);
							}
							result.push_back(new_centerline);
						}
					}
				}
			}
		}
		return result;
	}

	/// Split the fiber at the specified index. If the index is an end point, only one fiber is returned
	std::vector< stim::centerline<T> > split(unsigned int idx){

		std::vector< stim::centerline<T> > fl;				//create an array to store up to two fibers
		size_t N = size();

		//if the index is an end point, only the existing fiber is returned
		if(idx == 0 || idx == N-1){
			fl.resize(1);							//set the size of the fiber to 1
			fl[0] = *this;							//copy the current fiber
		}

		//if the index is not an end point
		else{

			unsigned int N1 = idx + 1;					//calculate the size of both fibers
			unsigned int N2 = N - idx;

			fl.resize(2);								//set the array size to 2

			fl[0] = stim::centerline<T>(N1);			//set the size of each fiber
			fl[1] = stim::centerline<T>(N2);

			//copy both halves of the fiber
			unsigned int i;

			//first half
			for(i = 0; i < N1; i++)					//for each centerline point
				fl[0][i] = std::vector< stim::vec3<T> >::at(i);
			fl[0].init();							//initialize the length vector

			//second half
			for(i = 0; i < N2; i++)
				fl[1][i] = std::vector< stim::vec3<T> >::at(idx+i);
			fl[1].init();							//initialize the length vector
		}

		return fl;										//return the array

	}

	/// Outputs the fiber as a string
	std::string str(){
		std::stringstream ss;
		size_t N = std::vector< stim::vec3<T> >::size();
		ss << "---------[" << N << "]---------" << std::endl;
		for (size_t i = 0; i < N; i++)
			ss << std::vector< stim::vec3<T> >::at(i) << std::endl;
		ss << "--------------------" << std::endl;

		return ss.str();
	}

	/// Back method returns the last point in the fiber
	stim::vec3<T> back(){
		return std::vector< stim::vec3<T> >::back();
	}

		////resample a fiber in the network
	stim::centerline<T> resample(T spacing)
	{	
		//std::cout<<"fiber::resample()"<<std::endl;

		stim::vec3<T> v;    //v-direction vector of the segment
		stim::vec3<T> p;      //- intermediate point to be added
		stim::vec3<T> p1;   // p1 - starting point of an segment on the fiber,
		stim::vec3<T> p2;   // p2 - ending point,
		//double sum=0;  //distance summation

		size_t N = size();

		centerline<T> new_c; // initialize list of new resampled points on the fiber
		// for each point on the centerline (skip if it is the last point on centerline)
		for(unsigned int f=0; f< N-1; f++)
		{			
			p1 = at(f); 
			p2 = at(f+1);
			v = p2 - p1;
			
			T lengthSegment = v.len();			//find Length of the segment as distance between the starting and ending points of the segment

			if(lengthSegment >= spacing){ // if length of the segment is greater than standard deviation resample
				
				// repeat resampling until accumulated stepsize is equsl to length of the segment
				for(T step=0.0; step<lengthSegment; step+=spacing){
					// calculate the resampled point by travelling step size in the direction of normalized gradient vector
					p = p1 + v * (step / lengthSegment);
					
					// add this resampled points to the new fiber list
					new_c.push_back(p);
				}
			}
			else       // length of the segment is now less than standard deviation, push the ending point of the segment and proceed to the next point in the fiber
				new_c.push_back(at(f));
		}
		new_c.push_back(at(N-1));   //add the last point on the fiber to the new fiber list
		//centerline newFiber(newPointList);
		return new_c;
	}

};



}	//end namespace stim



#endif
