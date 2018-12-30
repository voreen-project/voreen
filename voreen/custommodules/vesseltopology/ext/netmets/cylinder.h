#ifndef STIM_CYLINDER_H
#define STIM_CYLINDER_H
#include <iostream>
#include <stim/math/circle.h>
#include "centerline.h"
#include <stim/visualization/obj.h>


namespace stim
{
template<typename T>
class cylinder : public centerline<T> {
protected:
	
	using stim::centerline<T>::d;

	std::vector< stim::vec3<T> > U;					//stores the array of U vectors defining the Frenet frame
	std::vector< T > R;				//stores a list of magnitudes for each point in the centerline (assuming mags[0] is the radius)

	using stim::centerline<T>::findIdx;

	void init() {
		U.resize(size());			//allocate space for the frenet frame vectors
//		if (R.size() != 0)
		R.resize(size());

		stim::circle<T> c;								//create a circle
		stim::vec3<T> d0, d1;
		for (size_t i = 0; i < size() - 1; i++) {		//for each line segment in the centerline
			c.rotate(d(i));								//rotate the circle to match that normal
			U[i] = c.U;									//save the U vector from the circle
		}
		U[size() - 1] = c.U;							//for the last point, duplicate the final frenet frame vector
	}
	
	//calculates the U values for each point to initialize the frenet frame
	//	this function assumes that the centerline has already been set

public:

	using stim::centerline<T>::size;
	using stim::centerline<T>::at;
	using stim::centerline<T>::L;
	using stim::centerline<T>::length;

	cylinder() : centerline<T>(){}

	cylinder(centerline<T>c) : centerline<T>(c) {
		init();
	}

	cylinder(std::vector<stim::vec3<T> > p, std::vector<T> s)
		: centerline<T>(p)
	{
		R = s;
		init();
	}

	cylinder(stim::centerline<T> p, std::vector<T> s)
	{
		//was d = s;
		p = s;
		init();
	}
	
	//cylinder(centerline<T>c, T r) : centerline(c) {
	//	init();
	//	//add_mag(r);
	//}

	//initialize a cylinder with a list of points and magnitude values
	//cylinder(centerline<T>c, std::vector<T> r) : centerline(c){
	//	init();
	//	//add_mag(r);
	//}

	//copy the original radius
	void copy_r(std::vector<T> radius) {
		for (unsigned i = 0; i < radius.size(); i++)
			R[i] = radius[i];
	}

	///Returns magnitude i at the given length into the fiber (based on the pvalue).
	///Interpolates the position along the line.
	///@param l: the location of the in the cylinder.
	///@param idx: integer location of the point closest to l but prior to it.
	T r(T l, int idx) {
		T a = (l - L[idx]) / (L[idx + 1] - L[idx]);
		T v2 = (R[idx] + (R[idx + 1] - R[idx])*a);
		
		return v2;
	}

	///Returns the ith magnitude at the given p-value (p value ranges from 0 to 1).
	///interpolates the radius along the line.
	///@param pvalue: the location of the in the cylinder, from 0 (beginning to 1).
	T rl(T pvalue) {
		if (pvalue <= 0.0) return R[0];
		if (pvalue >= 1.0) return R[size() - 1];

		T l = pvalue*L[L.size() - 1];
		int idx = findIdx(l);
		return r(l, idx);
	}

	///	Returns the magnitude at the given index
	///	@param i is the index of the desired point
	/// @param r is the index of the magnitude value
	T r(unsigned i) {
		return R[i];
	}


	///adds a magnitude to each point in the cylinder
	/*void add_mag(V val = 0) {
		if (M.size() == 0) M.resize(size());	//if the magnitude vector isn't initialized, resize it to match the centerline
		for (size_t i = 0; i < size(); i++)		//for each point
			R[i].push_back(val);				//add this value to the magnitude vector at each point
	}*/

	//adds a magnitude based on a list of magnitudes for each point
	/*void add_mag(std::vector<T> val) {
		if (M.size() == 0) M.resize(size());	//if the magnitude vector isn't initialized, resize it to match the centerline
		for (size_t i = 0; i < size(); i++)		//for each point
			R[i].push_back(val[i]);				//add this value to the magnitude vector at each point
	}*/

	//sets the value of magnitude m at point i
	void set_r(size_t i, T r) {
		R[i] = r;
	}

	/*size_t nmags() {
		if (M.size() == 0) return 0;
		else return R[0].size();
	}*/

	///Returns a circle representing the cylinder cross section at point i
	stim::circle<T> circ(size_t i) {
		return stim::circle<T>(at(i), R[i], d(i), U[i]);
	}

	///Returns an OBJ object representing the cylinder with a radial tesselation value of N using magnitude m
	stim::obj<T> OBJ(size_t N) {
		stim::obj<T> out;								//create an OBJ object
		stim::circle<T> c0, c1;
		std::vector< stim::vec3<T> > p0, p1;			//allocate space for the point sets representing the circles bounding each cylinder segment
		T u0, u1, v0, v1;											//allocate variables to store running texture coordinates
		for (size_t i = 1; i < size(); i++) {			//for each line segment in the cylinder
			c0 = circ(i - 1);							//get the two circles bounding the line segment
			c1 = circ(i);

			p0 = c0.points(N);							//get t points for each of the end caps
			p1 = c1.points(N);

			u0 = L[i - 1] / length();						//calculate the texture coordinate (u, v) where u runs along the cylinder length
			u1 = L[i] / length();
				
			for (size_t n = 1; n < N; n++) {				//for each point in the circle
				v0 = (double)(n-1) / (double)(N - 1);			//v texture coordinate runs around the cylinder
				v1 = (double)(n) / (double)(N - 1);
				out.Begin(OBJ_FACE);						//start a face (quad)
					out.TexCoord(u0, v0);
					out.Vertex(p0[n - 1][0], p0[n - 1][1], p0[n - 1][2]);	//output the points composing a strip of quads wrapping the cylinder segment
					out.TexCoord(u1, v0);
					out.Vertex(p1[n - 1][0], p1[n - 1][1], p1[n - 1][2]);
				
					out.TexCoord(u0, v1);
					out.Vertex(p1[n][0], p1[n][1], p1[n][2]);
					out.TexCoord(u1, v1);
					out.Vertex(p0[n][0], p0[n][1], p0[n][2]);
				out.End();
			}
			v0 = (double)(N - 2) / (double)(N - 1);			//v texture coordinate runs around the cylinder
			v1 = 1.0;
			out.Begin(OBJ_FACE);
				out.TexCoord(u0, v0);
				out.Vertex(p0[N - 1][0], p0[N - 1][1], p0[N - 1][2]);	//output the points composing a strip of quads wrapping the cylinder segment
				out.TexCoord(u1, v0);
				out.Vertex(p1[N - 1][0], p1[N - 1][1], p1[N - 1][2]);

				out.TexCoord(u0, v1);
				out.Vertex(p1[0][0], p1[0][1], p1[0][2]);
				out.TexCoord(u1, v1);
				out.Vertex(p0[0][0], p0[0][1], p0[0][2]);
			out.End();
		}
		return out;
	}

	std::string str() {
		std::stringstream ss;
		size_t N = std::vector< stim::vec3<T> >::size();
		ss << "---------[" << N << "]---------" << std::endl;
		for (size_t i = 0; i < N; i++)
			ss << std::vector< stim::vec3<T> >::at(i) << "   r = " << R[i] << "   u = " << U[i] << std::endl;
		ss << "--------------------" << std::endl;

		return ss.str();
	}

	/// Integrates a magnitude value along the cylinder.
	/// @param m is the magnitude value to be integrated (this is usually the radius)
	T integrate() {
		T sum = 0;							//initialize the integral to zero
		if (L.size() == 1)
			return sum;
		else {					
			T m0, m1;						//allocate space for both magnitudes in a single segment
			m0 = R[0];					//initialize the first point and magnitude to the first point in the cylinder
			T len = L[1];					//allocate space for the segment length


			for (unsigned p = 1; p < size(); p++) {				//for every consecutive point in the cylinder
				m1 = R[p];
				if (p > 1) len = (L[p] - L[p - 1]);		//calculate the segment length using the L array
				sum += (m0 + m1) / (T)2.0 * len;				//add the average magnitude, weighted by the segment length
				m0 = m1;									//move to the next segment by shifting points
			}
			return sum;			//return the integral
		}
	}

	/// Resamples the cylinder to provide a maximum distance of "spacing" between centerline points. All current
	///		centerline points are guaranteed to exist in the new cylinder
	/// @param spacing is the maximum spacing allowed between sample points
	cylinder<T> resample(T spacing) {
		cylinder<T> c = stim::centerline<T>::resample(spacing);			//resample the centerline and use it to create a new cylinder

		//size_t nm = nmags();											//get the number of magnitude values in the current cylinder
		//if (nm > 0) {													//if there are magnitude values
		//	std::vector<T> magvec(nm, 0);							//create a magnitude vector for a single point
		//	c.M.resize(c.size(), magvec);									//allocate space for a magnitude vector at each point of the new cylinder
		//}

		T l, t;
		for (size_t i = 0; i < c.size(); i++) {							//for each point in the new cylinder
			l = c.L[i];													//get the length along the new cylinder
			t = l / length();										//calculate the parameter value along the new cylinder
			//for (size_t mag = 0; mag < nm; mag++) {							//for each magnitude value
			c.R[i] = r(t);								//retrieve the interpolated magnitude from the current cylinder and store it in the new one
			//}
		}
		return c;
	}

	std::vector< stim::cylinder<T> > split(unsigned int idx) {

		unsigned N = size();
		std::vector< stim::centerline<T> > LL;
		LL.resize(2);
		LL = (*this).centerline<T>::split(idx);
		std::vector< stim::cylinder<T> > C(LL.size());
		unsigned i = 0;
		C[0] = LL[0];
		//C[0].R.resize(idx);
		for (; i < idx + 1; i++) {
			//for(unsigned d = 0; d < 3; d++)
			//C[0][i][d] = LL[0].c[i][d];
			C[0].R[i] = R[i];
			//C[0].R[i].resize(1);
		}
		if (C.size() == 2) {
			C[1] = LL[1];
			i--;
			//C[1].M.resize(N - idx);
			for (; i < N; i++) {
				//for(unsigned d = 0; d < 3; d++)
				//C[1][i - idx][d] = LL[1].c[i - idx][d];
				C[1].R[i - idx] = R[i];
				//C[1].M[i - idx].resize(1);
			}
		}

		return C;
	}


		/*
		///inits the cylinder from a list of points (std::vector of stim::vec3 --inP) and magnitudes (inM)
		void
		init(centerline inP, std::vector< std::vector<T> > inM) {
			M = inM;									//the magnitude vector can be copied directly
			(*this) = inP;								//the centerline can be copied to this class directly
			stim::vec3<float> v1;
			stim::vec3<float> v2;
			e.resize(inP.size());

			norms.resize(inP.size());
			Us.resize(inP.size());

			if(inP.size() < 2)
				return;

			//calculate each L.
			L.resize(inP.size());						//the number of precomputed lengths will equal the number of points
			T temp = (T)0;								//length up to that point
			L[0] = temp;
			for(size_t i = 1; i < L.size(); i++)
			{
				temp += (inP[i-1] - inP[i]).len();
				L[i] = temp;
			}

			stim::vec3<T> dr = (inP[1] - inP[0]).norm();
			s = stim::circle<T>(inP[0], inR[0][0], dr, stim::vec3<T>(1,0,0));
			norms[0] = s.N;
			e[0] = s;
			Us[0] = e[0].U;
			for(size_t i = 1; i < inP.size()-1; i++)
			{
				s.center(inP[i]);
				v1 = (inP[i] - inP[i-1]).norm();
				v2 = (inP[i+1] - inP[i]).norm();
				dr = (v1+v2).norm();
				s.normal(dr);

				norms[i] = s.N;

				s.scale(inR[i][0]/inR[i-1][0]);
				e[i] = s;
				Us[i] = e[i].U;
			}
			
			int j = inP.size()-1;
			s.center(inP[j]);
			dr = (inP[j] - inP[j-1]).norm();
			s.normal(dr);

			norms[j] = s.N;

			s.scale(inR[j][0]/inR[j-1][0]);
			e[j] = s;
			Us[j] = e[j].U;
		}
		
		///returns the direction vector at point idx.
		stim::vec3<T>
		d(int idx)
		{
			if(idx == 0)
			{
				stim::vec3<T> temp(
						c[idx+1][0]-c[idx][0],	
						c[idx+1][1]-c[idx][1],	
						c[idx+1][2]-c[idx][2]
						);	
//				return (e[idx+1].P - e[idx].P).norm();
				return (temp.norm());
			}
			else if(idx == N-1)
			{
				stim::vec3<T> temp(
						c[idx][0]-c[idx+1][0],	
						c[idx][1]-c[idx+1][1],	
						c[idx][2]-c[idx+1][2]
						);	
			//	return (e[idx].P - e[idx-1].P).norm();
				return (temp.norm());
			}
			else
			{
//				return (e[idx+1].P - e[idx].P).norm();
//				stim::vec3<float> v1 = (e[idx].P-e[idx-1].P).norm();
				stim::vec3<T> v1(
						c[idx][0]-c[idx-1][0],	
						c[idx][1]-c[idx-1][1],	
						c[idx][2]-c[idx-1][2]
						);	
						
//				stim::vec3<float> v2 = (e[idx+1].P-e[idx].P).norm();
				stim::vec3<T> v2(
						c[idx+1][0]-c[idx][0],	
						c[idx+1][1]-c[idx][1],	
						c[idx+1][2]-c[idx][2]
						);
					
				return (v1.norm()+v2.norm()).norm();			
			} 
	//		return e[idx].N;	

		}

		stim::vec3<T>
		d(T l, int idx)
		{
			if(idx == 0 || idx == N-1)
			{
				return norms[idx];
//				return e[idx].N;
			}
			else
			{
				
				T rat = (l-L[idx])/(L[idx+1]-L[idx]);
				return(	norms[idx] + (norms[idx+1] - norms[idx])*rat);
//				return(	e[idx].N + (e[idx+1].N - e[idx].N)*rat);
			} 
		}


		///finds the index of the point closest to the length l on the lower bound.
		///binary search.
		int
		findIdx(T l)
		{
			unsigned int i = L.size()/2;
			unsigned int max = L.size()-1;
			unsigned int min = 0;
			while(i > 0 && i < L.size()-1)
			{
//				std::cerr << "Trying " << i << std::endl;
//				std::cerr << "l is " << l << ", L[" << i << "]" << L[i] << std::endl;
				if(l < L[i])
				{
					max = i;
					i = min+(max-min)/2;
				}
				else if(L[i] <= l && L[i+1] >= l)
				{
					break;
				}
				else
				{
					min = i;
					i = min+(max-min)/2;
				}
			}
			return i;
		}

	public:
		///default constructor
		cylinder()
		// : centerline<T>()
		{

		}

		///constructor to create a cylinder from a set of points, radii, and the number of sides for the cylinder.
		///@param inP:  Vector of stim vec3 composing the points of the centerline.
		///@param inM:  Vector of stim vecs composing the radii of the centerline.
		cylinder(std::vector<stim::vec3<T> > inP, std::vector<stim::vec<T> > inM)
			: centerline<T>(inP)
		{
			init(inP, inM);
		}

		///constructor to create a cylinder from a set of points, radii, and the number of sides for the cylinder.
		///@param inP:  Vector of stim vec3 composing the points of the centerline.
		///@param inM:  Vector of stim vecs composing the radii of the centerline.
		cylinder(std::vector<stim::vec3<T> > inP, std::vector< T > inM)
			: centerline<T>(inP)
		{
			std::vector<stim::vec<T> > temp;
			stim::vec<T> zero(0.0,0.0);
			temp.resize(inM.size(), zero);
			for(int i = 0; i < inM.size(); i++)
				temp[i][0] = inR[i];
			init(inP, temp);
		}


		///Constructor defines a cylinder with centerline inP and magnitudes of zero
		///@param inP: Vector of stim vec3 composing the points of the centerline
		cylinder(std::vector< stim::vec3<T> > inP)
			: centerline<T>(inP)
		{
			std::vector< stim::vec<T> > inM;						//create an array of arbitrary magnitudes

			stim::vec<T> zero;
			zero.push_back(0);

			inM.resize(inP.size(), zero);								//initialize the magnitude values to zero
			init(inP, inM);
		}

		//assignment operator creates a cylinder from a centerline (default radius is zero)
		cylinder& operator=(stim::centerline<T> c) {
			init(c);
		}

		///Returns the number of points on the cylinder centerline

		unsigned int size(){
			return N;
		}

		
		///Returns a position vector at the given p-value (p value ranges from 0 to 1).
		///interpolates the position along the line.
		///@param pvalue: the location of the in the cylinder, from 0 (beginning to 1).
		stim::vec3<T>
		p(T pvalue)
		{
			if(pvalue < 0.0 || pvalue > 1.0)
			{
				return stim::vec3<float>(-1,-1,-1);
			}
			T l = pvalue*L[L.size()-1];
			int idx = findIdx(l);
			return (p(l,idx));
		}

		///Returns a position vector at the given length into the fiber (based on the pvalue).
		///Interpolates the radius along the line.
		///@param l: the location of the in the cylinder.
		///@param idx: integer location of the point closest to l but prior to it.
		stim::vec3<T>
		p(T l, int idx)
		{
				T rat = (l-L[idx])/(L[idx+1]-L[idx]);
				stim::vec3<T> v1(
						c[idx][0],	
						c[idx][1],	
						c[idx][2]
						);	
						
				stim::vec3<T> v2(
						c[idx+1][0],	
						c[idx+1][1],	
						c[idx+1][2]
						);
//			return(	e[idx].P + (e[idx+1].P-e[idx].P)*rat);
			return(	v1 + (v2-v1)*rat);
//			return(
//			return (pos[idx] + (pos[idx+1]-pos[idx])*((l-L[idx])/(L[idx+1]- L[idx])));
		}

		///Returns a radius at the given p-value (p value ranges from 0 to 1).
		///interpolates the radius along the line.
		///@param pvalue: the location of the in the cylinder, from 0 (beginning to 1).
		T
		r(T pvalue)
		{
			if(pvalue < 0.0 || pvalue > 1.0){
				std::cerr<<"Error, value "<<pvalue<<" is outside of [0 1]."<<std::endl;
				exit(1);
			}
			T l = pvalue*L[L.size()-1];
			int idx = findIdx(l);
			return (r(l,idx));
		}

		///Returns a radius at the given length into the fiber (based on the pvalue).
		///Interpolates the position along the line.
		///@param l: the location of the in the cylinder.
		///@param idx: integer location of the point closest to l but prior to it.
		T
		r(T l, int idx)
		{
				T rat = (l-L[idx])/(L[idx+1]-L[idx]);
			T v1 = (e[idx].U.len() + (e[idx+1].U.len() - e[idx].U.len())*rat);
			T v3 = (Us[idx].len() + (Us[idx+1].len() - Us[idx].len())*rat);
			T v2 = (mags[idx][0] + (mags[idx+1][0]-mags[idx][0])*rat);
//			std::cout << (float)v1 = (float) v2 << std::endl;
			if(abs(v3 - v1) >= 10e-6)
			{
				std::cout << "-------------------------" << std::endl;
				std::cout << e[idx].str() << std::endl << std::endl;
				std::cout << Us[idx].str() << std::endl;
				std::cout << (float)v1 - (float) v2 << std::endl;
				std::cout << "failed" << std::endl;
			}
//			std::cout << e[idx].U.len() << " " << mags[idx][0] << std::endl;
//			std::cout << v2 << std::endl;
			return(v2);
//			return (mags[idx][0] + (mags[idx+1][0]-mags[idx][0])*rat);
	//	(
		}

		///	Returns the magnitude at the given index
		///	@param i is the index of the desired point
		/// @param m is the index of the magnitude value
		T ri(unsigned i, unsigned m = 0){
			return mags[i][m];
		}

		/// Adds a new magnitude value to all points
		/// @param m is the starting value for the new magnitude
		void add_mag(T m = 0){
			for(unsigned int p = 0; p < N; p++)
				mags[p].push_back(m);
		}

		/// Sets a magnitude value
		/// @param val is the new value for the magnitude
		/// @param p is the point index for the magnitude to be set
		/// @param m is the index for the magnitude
		void set_mag(T val, unsigned p, unsigned m = 0){
			mags[p][m] = val;
		}

		/// Returns the number of magnitude values at each point
		unsigned nmags(){
			return mags[0].size();
		}

		///returns the position of the point with a given pvalue and theta on the surface
		///in x, y, z coordinates. Theta is in degrees from 0 to 360.
		///@param pvalue: the location of the in the cylinder, from 0 (beginning to 1).
		///@param theta: the angle to the point of a circle.
		stim::vec3<T>
		surf(T pvalue, T theta)
		{
			if(pvalue < 0.0 || pvalue > 1.0)
			{
				return stim::vec3<float>(-1,-1,-1);
			} else {
			T l = pvalue*L[L.size()-1];
			int idx = findIdx(l);
			stim::vec3<T> ps = p(l, idx); 
			T m = r(l, idx);
			s = e[idx];
			s.center(ps);
			s.normal(d(l, idx));
			s.scale(m/e[idx].U.len());
			return(s.p(theta));
			}
		}

		///returns a vector of points necessary to create a circle at every position in the fiber.
		///@param sides: the number of sides of each circle.	
		std::vector<std::vector<vec3<T> > >
		getPoints(int sides)
		{
			std::vector<std::vector <vec3<T> > > points;
			points.resize(N);
			for(int i = 0; i < N; i++)
			{
				points[i] = e[i].getPoints(sides);
			}
			return points;
		}

		///returns the total length of the line at index j.
		T
		getl(int j)
		{
			return (L[j]);
		}
		/// Allows a point on the centerline to be accessed using bracket notation

		vec3<T> operator[](unsigned int i){
			return e[i].P;
		}

		/// Returns the total length of the cylinder centerline
		T length(){
			return L.back();
		}

		/// Integrates a magnitude value along the cylinder.
		/// @param m is the magnitude value to be integrated (this is usually the radius)
		T integrate(unsigned m = 0){

			T M = 0;						//initialize the integral to zero
			T m0, m1;						//allocate space for both magnitudes in a single segment

			//vec3<T> p0, p1;					//allocate space for both points in a single segment

			m0 = mags[0][m];				//initialize the first point and magnitude to the first point in the cylinder
			//p0 = pos[0];

			T len = L[0];						//allocate space for the segment length

			//for every consecutive point in the cylinder
			for(unsigned p = 1; p < N; p++){

				//p1 = pos[p];							//get the position and magnitude for the next point
				m1 = mags[p][m];

				if(p > 1) len = (L[p-1] - L[p-2]);		//calculate the segment length using the L array

				//add the average magnitude, weighted by the segment length
				M += (m0 + m1)/(T)2.0 * len;

				m0 = m1;								//move to the next segment by shifting points
			}
			return M;			//return the integral
		}

		/// Averages a magnitude value across the cylinder
		/// @param m is the magnitude value to be averaged (this is usually the radius)
		T average(unsigned m = 0){			

			//return the average magnitude
			return integrate(m) / L.back();
		}

		/// Resamples the cylinder to provide a maximum distance of "spacing" between centerline points. All current
		///		centerline points are guaranteed to exist in the new cylinder
		/// @param spacing is the maximum spacing allowed between sample points
		cylinder<T> resample(T spacing){

			std::vector< vec3<T> > result;

			vec3<T> p0 = e[0].P;									//initialize p0 to the first point on the centerline
			vec3<T> p1;
			unsigned N = size();									//number of points in the current centerline

			//for each line segment on the centerline
			for(unsigned int i = 1; i < N; i++){
				p1 = e[i].P;										//get the second point in the line segment

				vec3<T> v = p1 - p0;								//calculate the vector between these two points
				T d = v.len();										//calculate the distance between these two points (length of the line segment)

				size_t nsteps = (size_t)std::ceil(d / spacing);		//calculate the number of steps to take along the segment to meet the spacing criteria
				T stepsize = (T)1.0 / nsteps;						//calculate the parametric step size between new centerline points

				//for each step along the line segment
				for(unsigned s = 0; s < nsteps; s++){
					T alpha = stepsize * s;							//calculate the fraction of the distance along the line segment covered
					result.push_back(p0 + alpha * v);				//push the point at alpha position along the line segment
				}

				p0 = p1;											//shift the points to move to the next line segment
			}

			result.push_back(e[size() - 1].P);						//push the last point in the centerline

			return cylinder<T>(result);

		}*/

		
};

}
#endif
