#ifndef STIM_GL_NETWORK
#define STIM_GL_NETWORK

#include <GL/glut.h>
#include "network.h"
#include <stim/visualization/aaboundingbox.h>

namespace stim{

template <typename T>
class gl_network : public stim::network<T>{

protected:
	using stim::network<T>::E;
	using stim::network<T>::V;

	GLuint dlist;

public:

	/// Default constructor
	gl_network() : stim::network<T>(){
		dlist = 0;
	}

	/// Constructor creates a gl_network from a stim::network
	gl_network(stim::network<T> N) : stim::network<T>(N){
		dlist = 0;
	}

	/// Fills the parameters with the minimum and maximum spatial positions in the network,
	///     specifying a bounding box for the network geometry
	aaboundingbox<T> boundingbox(){

		aaboundingbox<T> bb;								//create a bounding box

		//loop through every edge
		for(unsigned e = 0; e < E.size(); e++){
			//loop through every point
			for(unsigned p = 0; p < E[e].size(); p++)
				bb.expand(E[e][p]);						//expand the bounding box to include the point
		}

		return bb;								//return the bounding box
	}

	///render cylinder based on points from the top/bottom hat
	///@param C1 set of points from one of the hat
	void renderCylinder(std::vector< stim::vec3<T> > C1, std::vector< stim::vec3<T> > C2, stim::vec3<T> P1, stim::vec3<T> P2) {
		glBegin(GL_QUAD_STRIP);
		for (unsigned i = 0; i < C1.size(); i++) {				// for every point on the circle
			stim::vec3<T> n1 = C1[i] - P1;						// set normal vector for every vertex on the quads
			stim::vec3<T> n2 = C2[i] - P2;
			n1 = n1.norm();
			n2 = n2.norm();
			glNormal3f(n1[0], n1[1], n1[2]);
			glVertex3f(C1[i][0], C1[i][1], C1[i][2]);
			glNormal3f(n2[0], n2[1], n2[2]);
			glVertex3f(C2[i][0], C2[i][1], C2[i][2]);											
		}	
		glEnd();
	}

	///render the vertex as sphere based on glut build-in function
	///@param x, y, z are the three coordinates of the center point
	///@param radius is the radius of the sphere
	///@param subdivisions is the slice/stride along/around z-axis
	void renderBall(T x, T y, T z, T radius, int subdivisions) {
		glPushMatrix();
		glTranslatef(x, y, z);
		glutSolidSphere(radius, subdivisions, subdivisions);
		glPopMatrix();
	}

	///render the vertex as sphere based on transformation
	///@param x, y, z are the three coordinates of the center point
	///@param radius is the radius of the sphere
	///@param slice is the number of subdivisions around the z-axis
	///@param stack is the number of subdivisions along the z-axis
	void renderBall(T x, T y, T z, T radius, T slice, T stack) {
		T step_z = stim::PI / slice;					// step angle along z-axis
		T step_xy = 2 * stim::PI / stack;				// step angle in xy-plane
		T xx[4], yy[4], zz[4];							// store coordinates

		T angle_z = 0.0;								// start angle
		T angle_xy = 0.0;

		glBegin(GL_QUADS);
		for (unsigned i = 0; i < slice; i++) {			// around the z-axis
			angle_z = i * step_z;						// step step_z each time

			for (unsigned j = 0; j < stack; j++) {		// along the z-axis
				angle_xy = j * step_xy;					// step step_xy each time, draw floor by floor

				xx[0] = radius * std::sin(angle_z) * std::cos(angle_xy);	// four vertices
				yy[0] = radius * std::sin(angle_z) * std::sin(angle_xy);
				zz[0] = radius * std::cos(angle_z);

				xx[1] = radius * std::sin(angle_z + step_z) * std::cos(angle_xy);
				yy[1] = radius * std::sin(angle_z + step_z) * std::sin(angle_xy);
				zz[1] = radius * std::cos(angle_z + step_z);

				xx[2] = radius * std::sin(angle_z + step_z) * std::cos(angle_xy + step_xy);
				yy[2] = radius * std::sin(angle_z + step_z) * std::sin(angle_xy + step_xy);
				zz[2] = radius * std::cos(angle_z + step_z);

				xx[3] = radius * std::sin(angle_z) * std::cos(angle_xy + step_xy);
				yy[3] = radius * std::sin(angle_z) * std::sin(angle_xy + step_xy);
				zz[3] = radius * std::cos(angle_z);

				for (unsigned k = 0; k < 4; k++) {
					glVertex3f(x + xx[k], y + yy[k], z + zz[k]);			// draw the floor plane
				}
			}
		}
		glEnd();
	}

	/// Render the network centerline as a series of line strips.
	/// glCenterline0 is for only one input
	void glCenterline0(){
		if (!glIsList(dlist)) {					//if dlist isn't a display list, create it
			dlist = glGenLists(1);				//generate a display list
			glNewList(dlist, GL_COMPILE);		//start a new display list
			for (unsigned e = 0; e < E.size(); e++) {				//for each edge in the network
				glBegin(GL_LINE_STRIP);
				for (unsigned p = 0; p < E[e].size(); p++) {			//for each point on that edge
					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);		//set the vertex position based on the current point
					glTexCoord1f(0);									//set white color
				}
				glEnd();
			}
			glEndList();						//end the display list
		}
		glCallList(dlist);					// render the display list
	}

	///render the network centerline as a series of line strips(when loading at least two networks, otherwise using glCenterline0())
	///colors are based on metric values
	void glCenterline(){

		if(!glIsList(dlist)){					//if dlist isn't a display list, create it
			dlist = glGenLists(1);				//generate a display list
			glNewList(dlist, GL_COMPILE);		//start a new display list
			for(unsigned e = 0; e < E.size(); e++){				//for each edge in the network
				//unsigned errormag_id = E[e].nmags() - 1;
				glBegin(GL_LINE_STRIP);
				for(unsigned p = 0; p < E[e].size(); p++){				//for each point on that edge
					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);		//set the vertex position based on the current point
					glTexCoord1f(E[e].r(p));							//set the texture coordinate based on the specified magnitude index
				}
				glEnd();
			}
			glEndList();						//end the display list
		}		
		glCallList(dlist);						//render the display list
	}

	///render the network cylinder as a series of tubes(when only one network loaded) 
	void glCylinder0(T scale = 1.0f, bool undo = false) {

		float r1, r2;
		if (undo == true)
			glDeleteLists(dlist, 1);								// delete display list 
		if (!glIsList(dlist)) {										// if dlist isn't a display list, create it
			dlist = glGenLists(1);									// generate a display list
			glNewList(dlist, GL_COMPILE);							// start a new display list
			for (unsigned e = 0; e < E.size(); e++) {				// for each edge in the network
				for (unsigned p = 1; p < E[e].size(); p++) {		// for each point on that edge
					stim::circle<T> C1 = E[e].circ(p - 1);
					stim::circle<T> C2 = E[e].circ(p);
					r1 = E[e].r(p - 1);
					r2 = E[e].r(p);
					C1.set_R(scale * r1);								// re-scale the circle to the same
					C2.set_R(scale * r2);
					std::vector< stim::vec3<T> > Cp1 = C1.glpoints(20);	// get 20 points on the circle plane
					std::vector< stim::vec3<T> > Cp2 = C2.glpoints(20);
					glBegin(GL_QUAD_STRIP);
					for (unsigned i = 0; i < Cp1.size(); i++) {
						glVertex3f(Cp1[i][0], Cp1[i][1], Cp1[i][2]);
						glVertex3f(Cp2[i][0], Cp2[i][1], Cp2[i][2]);
					}
					glEnd();
				}													// set the texture coordinate based on the specified magnitude index
			}
			for (unsigned n = 0; n < V.size(); n++) {
				for (unsigned i = 0; i < E.size(); i++) {
					if (E[i].v[0] == n) {
						r1 = E[i].r(0) * scale;
						break;
					}
					else if (E[i].v[1] == n) {
						r1 = E[i].r(E[i].size() - 1) * scale;
						break;
					}
				}
				renderBall(V[n][0], V[n][1], V[n][2], r1, 20);
			}
			glEndList();						// end the display list
		}
		glCallList(dlist);						// render the display list
	}

	///render the network cylinder as a series of tubes
	///colors are based on metric values 
	void glCylinder(float sigma, float radius) {

		if (radius != sigma)					// if render radius was changed by user, create a new display list
			glDeleteLists(dlist, 1);
		if (!glIsList(dlist)) {										// if dlist isn't a display list, create it
			dlist = glGenLists(1);									// generate a display list
			glNewList(dlist, GL_COMPILE);							// start a new display list
			for (unsigned e = 0; e < E.size(); e++) {				// for each edge in the network
				for (unsigned p = 1; p < E[e].size(); p++) {		// for each point on that edge
					stim::circle<T> C1 = E[e].circ(p - 1);
					stim::circle<T> C2 = E[e].circ(p);
					C1.set_R(2*radius);								// re-scale the circle to the same
					C2.set_R(2*radius);
					std::vector< stim::vec3<T> > Cp1 = C1.glpoints(20);// get 20 points on the circle plane
					std::vector< stim::vec3<T> > Cp2 = C2.glpoints(20);
					glBegin(GL_QUAD_STRIP);
					for (unsigned i = 0; i < Cp1.size(); i++) {
						stim::vec3<T> n1 = Cp1[i] - E[e][p - 1];	// set normal vector for every vertex on the quads
						stim::vec3<T> n2 = Cp2[i] - E[e][p];
						n1 = n1.norm();
						n2 = n2.norm();

						glNormal3f(n1[0], n1[1], n1[2]);
						glTexCoord1f(E[e].r(p - 1));
						glVertex3f(Cp1[i][0], Cp1[i][1], Cp1[i][2]);
						glNormal3f(n2[0], n2[1], n2[2]);
						glTexCoord1f(E[e].r(p));
						glVertex3f(Cp2[i][0], Cp2[i][1], Cp2[i][2]);
					}
					glEnd();
				}													// set the texture coordinate based on the specified magnitude index
			}
			for (unsigned n = 0; n < V.size(); n++) {
				size_t num = V[n].e[0].size();					// store the number of outgoing edge of that vertex
				if (num != 0) {									// if it has outgoing edge
					unsigned idx = V[n].e[0][0];				// find the index of first outgoing edge of that vertex
					glTexCoord1f(E[idx].r(0));					// bind the texture as metric of first point on that edge
				}					
				else {
					unsigned idx = V[n].e[1][0];				// find the index of first incoming edge of that vertex
					glTexCoord1f(E[idx].r(E[idx].size() - 1));	// bind the texture as metric of last point on that edge
				}
				renderBall(V[n][0], V[n][1], V[n][2], 2*radius, 20);
			}
			glEndList();						// end the display list
		}
		glCallList(dlist);						// render the display list
	}

	///render a T as a adjoint network of GT in transparancy(darkgreen, overlaid)
	void glAdjointCylinder(float sigma, float radius) {
		
		if (radius != sigma)					// if render radius was changed by user, create a new display list
			glDeleteLists(dlist + 4, 1);
		if (!glIsList(dlist + 4)) {
			glNewList(dlist + 4, GL_COMPILE);
			for (unsigned e = 0; e < E.size(); e++) {				// for each edge in the network
				for (unsigned p = 1; p < E[e].size(); p++) {		// for each point on that edge
					stim::circle<T> C1 = E[e].circ(p - 1);
					stim::circle<T> C2 = E[e].circ(p);
					C1.set_R(2 * radius);							// scale the circle to the same
					C2.set_R(2 * radius);
					std::vector< stim::vec3<T> >Cp1 = C1.glpoints(20);
					std::vector< stim::vec3<T> >Cp2 = C2.glpoints(20);
					glBegin(GL_QUAD_STRIP);
					for (unsigned i = 0; i < Cp1.size(); i++) {		// for every point on the circle(+1 means closing the circle)
						glVertex3f(Cp1[i][0], Cp1[i][1], Cp1[i][2]);
						glVertex3f(Cp2[i][0], Cp2[i][1], Cp2[i][2]);
					}
					glEnd();
				}								// set the texture coordinate based on the specified magnitude index
			}
			for (unsigned n = 0; n < V.size(); n++) {
				size_t num = V[n].e[0].size();					// store the number of outgoing edge of that vertex
				if (num != 0) {									// if it has outgoing edge
					unsigned idx = V[n].e[0][0];				// find the index of first outgoing edge of that vertex 
				}
				else {
					unsigned idx = V[n].e[1][0];				// find the index of first incoming edge of that vertex
				}
				renderBall(V[n][0], V[n][1], V[n][2], 2 * radius, 20);
			}
			glEndList();
		}
		glCallList(dlist + 4);
	}

	///render the network cylinder as series of tubes
	///@param I is a indicator: 0 -> GT, 1 -> T
	///@param map is the mapping relationship between two networks
	///@param colormap is the random generated color set for render
	void glRandColorCylinder(int I, std::vector<unsigned> map, std::vector<T> colormap, float sigma, float radius) {
		
		if (radius != sigma)									// if render radius was changed by user, create a new display list
			glDeleteLists(dlist + 2, 1);
		if (!glIsList(dlist + 2)) {								// if dlist isn't a display list, create it
			glNewList(dlist + 2, GL_COMPILE);					// start a new display list
			for (unsigned e = 0; e < E.size(); e++) {			// for each edge in the network
				if (map[e] != unsigned(-1)) {
					if (I == 0) {								// if it is to render GT
						glColor3f(colormap[e * 3 + 0], colormap[e * 3 + 1], colormap[e * 3 + 2]);
					}									
					else {										// if it is to render T
						glColor3f(colormap[map[e] * 3 + 0], colormap[map[e] * 3 + 1], colormap[map[e] * 3 + 2]);
					}

					for (unsigned p = 1; p < E[e].size(); p++) {// for each point on that edge
						stim::circle<T> C1 = E[e].circ(p - 1);
						stim::circle<T> C2 = E[e].circ(p);
						C1.set_R(2*radius);						// re-scale the circle to the same
						C2.set_R(2*radius);
						std::vector< stim::vec3<T> >Cp1 = C1.glpoints(20);
						std::vector< stim::vec3<T> >Cp2 = C2.glpoints(20);
						renderCylinder(Cp1, Cp2, E[e][p - 1], E[e][p]);
					}
				}
				else {
					glColor3f(1.0, 1.0, 1.0);					// white color for the un-mapping edges
					for (unsigned p = 1; p < E[e].size(); p++) {// for each point on that edge
						stim::circle<T> C1 = E[e].circ(p - 1);
						stim::circle<T> C2 = E[e].circ(p);
						C1.set_R(2*radius);						// scale the circle to the same
						C2.set_R(2*radius);
						std::vector< stim::vec3<T> >Cp1 = C1.glpoints(20);
						std::vector< stim::vec3<T> >Cp2 = C2.glpoints(20);
						renderCylinder(Cp1, Cp2, E[e][p - 1], E[e][p]);
					}
				}
			}
			for (unsigned n = 0; n < V.size(); n++) {
				size_t num_edge = V[n].e[0].size() + V[n].e[1].size();
				if (num_edge > 1) {					// if it is the joint vertex
					glColor3f(0.3, 0.3, 0.3);		// gray
					renderBall(V[n][0], V[n][1], V[n][2], 3*radius, 20);
				}
				else {								// if it is the terminal vertex
					glColor3f(0.6, 0.6, 0.6);		// more white gray
					renderBall(V[n][0], V[n][1], V[n][2], 3*radius, 20);
				}
			}
			glEndList();
		}
		glCallList(dlist + 2);
	}

	///render the network centerline as a series of line strips in random different color
	///@param I is a indicator: 0 -> GT, 1 -> T
	///@param map is the mapping relationship between two networks
	///@param colormap is the random generated color set for render
	void glRandColorCenterline(int I, std::vector<unsigned> map, std::vector<T> colormap) {
		if (!glIsList(dlist + 2)) {
			glNewList(dlist + 2, GL_COMPILE);
			for (unsigned e = 0; e < E.size(); e++) {
				if (map[e] != unsigned(-1)) {					// if it has corresponding edge in another network
					if (I == 0)									// if it is to render GT
						glColor3f(colormap[e * 3 + 0], colormap[e * 3 + 1], colormap[e * 3 + 2]);
					else										// if it is to render T
						glColor3f(colormap[map[e] * 3 + 0], colormap[map[e] * 3 + 1], colormap[map[e] * 3 + 2]);

					glBegin(GL_LINE_STRIP);
					for (unsigned p = 0; p < E[e].size(); p++) {
						glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
					}
					glEnd();
				}
				else {
					glColor3f(1.0, 1.0, 1.0);					// white color for the un-mapping edges
					glBegin(GL_LINE_STRIP);
					for (unsigned p = 0; p < E[e].size(); p++) {
						glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
					}
					glEnd();
				}
			}
			glEndList();
		}
		glCallList(dlist + 2);
	}

	void glAdjointCenterline() {
		if (!glIsList(dlist + 4)) {
			glNewList(dlist + 4, GL_COMPILE);
			for (unsigned e = 0; e < E.size(); e++) {				//for each edge in the network
																
				glBegin(GL_LINE_STRIP);
				for (unsigned p = 0; p < E[e].size(); p++) {			//for each point on that edge
					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);		//set the vertex position based on the current point
					glTexCoord1f(E[e].r(p));							//set the texture coordinate based on the specified magnitude index
				}
				glEnd();
			}
			glEndList();
		}
		glCallList(dlist + 4);
	}

	// highlight the difference part
	void glDifferenceCylinder(int I, std::vector<unsigned> map, std::vector<T> colormap, float sigma, float radius) {

		if (radius != sigma)									// if render radius was changed by user, create a new display list
			glDeleteLists(dlist + 6, 1);
		if (!glIsList(dlist + 6)) {								// if dlist isn't a display list, create it
			glNewList(dlist + 6, GL_COMPILE);					// start a new display list
			for (unsigned e = 0; e < E.size(); e++) {			// for each edge in the network
				if (map[e] != unsigned(-1)) {
					glEnable(GL_BLEND);									//enable color blend
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	//set blend function
					glDisable(GL_DEPTH_TEST);							//should disable depth

					if (I == 0) {								// if it is to render GT
						glColor4f(colormap[e * 3 + 0], colormap[e * 3 + 1], colormap[e * 3 + 2], 0.1);
					}
					else {										// if it is to render T
						glColor4f(colormap[map[e] * 3 + 0], colormap[map[e] * 3 + 1], colormap[map[e] * 3 + 2], 0.1);
					}

					for (unsigned p = 1; p < E[e].size(); p++) {// for each point on that edge
						stim::circle<T> C1 = E[e].circ(p - 1);
						stim::circle<T> C2 = E[e].circ(p);
						C1.set_R(2 * radius);						// re-scale the circle to the same
						C2.set_R(2 * radius);
						std::vector< stim::vec3<T> >Cp1 = C1.glpoints(20);
						std::vector< stim::vec3<T> >Cp2 = C2.glpoints(20);
						renderCylinder(Cp1, Cp2, E[e][p - 1], E[e][p]);
					}
					glDisable(GL_BLEND);
					glEnable(GL_DEPTH_TEST);
				}
				else {
					glColor3f(1.0, 1.0, 1.0);					// white color for the un-mapping edges
					for (unsigned p = 1; p < E[e].size(); p++) {// for each point on that edge
						stim::circle<T> C1 = E[e].circ(p - 1);
						stim::circle<T> C2 = E[e].circ(p);
						C1.set_R(2 * radius);						// scale the circle to the same
						C2.set_R(2 * radius);
						std::vector< stim::vec3<T> >Cp1 = C1.glpoints(20);
						std::vector< stim::vec3<T> >Cp2 = C2.glpoints(20);
						renderCylinder(Cp1, Cp2, E[e][p - 1], E[e][p]);
					}
				}
			}
			for (unsigned n = 0; n < V.size(); n++) {

				size_t num_edge = V[n].e[0].size() + V[n].e[1].size();
				if (num_edge > 1) {					// if it is the joint vertex
					glColor4f(0.3, 0.3, 0.3, 0.1);		// gray
					renderBall(V[n][0], V[n][1], V[n][2], 3 * radius, 20);
				}
				else {								// if it is the terminal vertex
					glColor4f(0.6, 0.6, 0.6, 0.1);		// more white gray
					renderBall(V[n][0], V[n][1], V[n][2], 3 * radius, 20);
				}
			}
			glEndList();
		}
		glCallList(dlist + 6);
	}

	//void glRandColorCenterlineGT(GLuint &dlist1, std::vector<unsigned> map, std::vector<T> colormap){
	//	if(!glIsList(dlist1)){
	//		dlist1 = glGenLists(1);
	//		glNewList(dlist1, GL_COMPILE);
	//		for(unsigned e = 0; e < E.size(); e++){
	//			if(map[e] != unsigned(-1)){
	//				glColor3f(colormap[e * 3 + 0], colormap[e * 3 + 1], colormap[e * 3 + 2]);
	//				glBegin(GL_LINE_STRIP);
	//				for(unsigned p = 0; p < E[e].size(); p++){
	//					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	//				}
	//				glEnd();
	//				for (unsigned p = 0; p < E[e].size() - 1; p++) {
	//					renderCylinder(E[e][p][0], E[e][p][1], E[e][p][2], E[e][p + 1][0], E[e][p + 1][1], E[e][p + 1][2], 10, 20);
	//				}
	//			}
	//			else{
	//				glColor3f(1.0, 1.0, 1.0);
	//				glBegin(GL_LINE_STRIP);
	//				for(unsigned p = 0; p < E[e].size(); p++){
	//					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	//				}
	//				glEnd();
	//			}
	//		}
	//		for (unsigned v = 0; v < V.size(); v++) {
	//			size_t num_edge = V[v].e[0].size() + V[v].e[1].size();
	//			if (num_edge > 1) {
	//				glColor3f(0.3, 0.3, 0.3);		// gray color for vertex
	//				renderBall(V[v][0], V[v][1], V[v][2], 20, 20);
	//			}
	//		}
	//		glEndList();
	//	}
	//	glCallList(dlist1);
	//}

	//void glRandColorCenterlineT(GLuint &dlist2, std::vector<unsigned> map, std::vector<T> colormap){
	//	if(!glIsList(dlist2)){
	//		dlist2 = glGenLists(1);
	//		glNewList(dlist2, GL_COMPILE);
	//		for(unsigned e = 0; e < E.size(); e++){
	//			if(map[e] != unsigned(-1)){
	//				glColor3f(colormap[map[e] * 3 + 0], colormap[map[e] * 3 + 1], colormap[map[e] * 3 + 2]);
	//				glBegin(GL_LINE_STRIP);
	//				for(unsigned p = 0; p < E[e].size(); p++){
	//					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	//				}
	//				glEnd();
	//				for (unsigned p = 0; p < E[e].size() - 1; p++) {
	//					renderCylinder(E[e][p][0], E[e][p][1], E[e][p][2], E[e][p + 1][0], E[e][p + 1][1], E[e][p + 1][2], 10, 20);
	//				}
	//			}
	//			else{
	//				glColor3f(1.0, 1.0, 1.0);
	//				glBegin(GL_LINE_STRIP);
	//				for(unsigned p = 0; p < E[e].size(); p++){
	//					glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	//				}
	//				glEnd();
	//			}
	//		}
	//		for (unsigned v = 0; v < V.size(); v++) {
	//			size_t num_edge = V[v].e[0].size() + V[v].e[1].size();
	//			if (num_edge > 1) {
	//				glColor3f(0.3, 0.3, 0.3);		// gray color for vertex
	//				renderBall(V[v][0], V[v][1], V[v][2], 20, 20);
	//			}
	//		}
	//		glEndList();
	//	}
	//	glCallList(dlist2);
	//}


	//void renderCylinder(T x1, T y1, T z1, T x2, T y2, T z2, T radius, int subdivisions) {
	//	T dx = x2 - x1;
	//	T dy = y2 - y1;
	//	T dz = z2 - z1;
	//	/// handle the degenerate case with an approximation
	//	if (dz == 0)
	//		dz = .00000001;
	//	T d = sqrt(dx*dx + dy*dy + dz*dz);					
	//	T ax = 57.2957795*acos(dz / d);						// 180°/pi
	//	if (dz < 0.0)
	//		ax = -ax;
	//	T rx = -dy*dz;
	//	T ry = dx*dz;

	//	glPushMatrix();
	//	glTranslatef(x1, y1, z1);
	//	glRotatef(ax, rx, ry, 0.0);

	//	glutSolidCylinder(radius, d, subdivisions, 1);
	//	glPopMatrix();
	//}


	/// render the network centerline from swc file as a series of strips in different colors based on the neuronal type
	/// glCenterline0_swc is for only one input
	/*void glCenterline0_swc() {
	if (!glIsList(dlist)) {						// if dlist isn't a display list, create it
	dlist = glGenLists(1);					// generate a display list
	glNewList(dlist, GL_COMPILE);			// start a new display list
	for (unsigned e = 0; e < E.size(); e++) {
	int type = NT[e];					// get the neuronal type
	switch (type) {
	case 0:
	glColor3f(1.0f, 1.0f, 1.0f);	// white for undefined
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 1:
	glColor3f(1.0f, 0.0f, 0.0f);	// red for soma
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 2:
	glColor3f(1.0f, 0.5f, 0.0f);	// orange for axon
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 3:
	glColor3f(1.0f, 1.0f, 0.0f);	// yellow for undefined
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 4:
	glColor3f(0.0f, 1.0f, 0.0f);	// green for undefined
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 5:
	glColor3f(0.0f, 1.0f, 1.0f);	// verdant for undefined
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 6:
	glColor3f(0.0f, 0.0f, 1.0f);	// blue for undefined
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	case 7:
	glColor3f(0.5f, 0.0f, 1.0f);	// purple for undefined
	glBegin(GL_LINE_STRIP);
	for (unsigned p = 0; p < E[e].size(); p++) {
	glVertex3f(E[e][p][0], E[e][p][1], E[e][p][2]);
	}
	glEnd();
	break;
	}
	}
	glEndList();						//end the display list
	}
	glCallList(dlist);					// render the display list
	}*/

	///render the network cylinder as a series of tubes
	///colors are based on metric values 
	//void glCylinder(float sigma) {
	//	if (!glIsList(dlist)) {					//if dlist isn't a display list, create it
	//		dlist = glGenLists(1);				//generate a display list
	//		glNewList(dlist, GL_COMPILE);		//start a new display list
	//		for (unsigned e = 0; e < E.size(); e++) {				//for each edge in the network
	//			for (unsigned p = 1; p < E[e].size() - 1; p++) {	// for each point on that edge
	//				stim::circle<T> C1 = E[e].circ(p - 1);
	//				stim::circle<T> C2 = E[e].circ(p);
	//				C1.set_R(2.5*sigma);							// scale the circle to the same
	//				C2.set_R(2.5*sigma);
	//				std::vector< stim::vec3<T> >Cp1 = C1.glpoints(20);
	//				std::vector< stim::vec3<T> >Cp2 = C2.glpoints(20);
	//				glBegin(GL_QUAD_STRIP);
	//				for (unsigned i = 0; i < Cp1.size(); i++) {		// for every point on the circle(+1 means closing the circle)
	//					glVertex3f(Cp1[i][0], Cp1[i][1], Cp1[i][2]);
	//					glVertex3f(Cp2[i][0], Cp2[i][1], Cp2[i][2]);
	//					glTexCoord1f(E[e].r(p));
	//				}
	//				glEnd();
	//			}								//set the texture coordinate based on the specified magnitude index
	//		}
	//		for (unsigned n = 0; n < V.size(); n++) {
	//			size_t num = V[n].e[0].size();					//store the number of outgoing edge of that vertex
	//			if (num != 0) {									//if it has outgoing edge
	//				unsigned idx = V[n].e[0][0];				//find the index of first outgoing edge of that vertex 
	//				glTexCoord1f(E[idx].r(0));					//bind the texture as metric of first point on that edge
	//			}
	//			else {
	//				unsigned idx = V[n].e[1][0];				//find the index of first incoming edge of that vertex
	//				glTexCoord1f(E[idx].r(E[idx].size() - 1));	//bind the texture as metric of last point on that edge
	//			}
	//			renderBall(V[n][0], V[n][1], V[n][2], 2.5*sigma, 20);
	//		}
	//		glEndList();						//end the display list
	//	}
	//	glCallList(dlist);						//render the display list
	//}

};		//end stim::gl_network class
};		//end stim namespace



#endif