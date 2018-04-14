/*! \file del_impl.cpp
        \brief The implementation of the 2D delaunay triangulation class.

        This class is a wrapper on the triangle package.
 */

#ifdef WIN32
#include <windows.h>
#endif

#include <iostream>
#include "../include/triangle_impl.hpp"
#include "../include/del_interface.hpp"
#include "../../../datastructures/plotrow.h"
#include <new>
#include <algorithm>
#include <utility>
#include <cmath>
#include <map>

#define REAL double

namespace tpp {

using std::cout;
using std::cerr;


/*** Voreen adaptations ***/

void Delaunay::generateVoronoiRegions(const voreen::PlotData& data, std::vector< std::list< tgt::dvec2 > >&  outVoronoiRegions) {
    //typedef struct triangulateio  TriangStruct;
    typedef struct triangulateio* pTriangStruct;

    pTriangStruct pvorout = (struct triangulateio*) vorout;
    
    // each input points leads to one voronoi region
    outVoronoiRegions.clear();
    outVoronoiRegions.resize(pvorout->edgesOfVoronoiFaces.size());
    
    for (size_t i=0; i<pvorout->edgesOfVoronoiFaces.size(); ++i) {
        double px = data.getRow(static_cast<int>(i)).getValueAt(0);
        double py = data.getRow(static_cast<int>(i)).getValueAt(1);
        
        // first use a map to efficiently sort vertices of voronoi regions by polar angle
        std::map<double, Point> vertexMap;
        
        // each voronoi region has some voronoi edges, we push the vertices of the edges into our sorted map
        for (std::vector<int>::const_iterator eit = pvorout->edgesOfVoronoiFaces[i].begin(); eit < pvorout->edgesOfVoronoiFaces[i].end(); ++eit) {
            int idx1 = pvorout->edgelist[*eit];
            double x1 = pvorout->pointlist[2*idx1];
            double y1 = pvorout->pointlist[2*idx1+1]; 
            
            int idx2 = pvorout->edgelist[(*eit)+1];
            double x2 = pvorout->pointlist[2*idx2];
            double y2 = pvorout->pointlist[2*idx2+1]; 
            
            // calculate polar angles
            double polar1 = atan2(y1-py, x1-px);
            double polar2 = atan2(y2-py, x2-px);
            
            // check wheter angle1 is already in map or not
            std::map<double, Point>::iterator lb = vertexMap.lower_bound(polar1);
            if (lb == vertexMap.end() || polar1 != lb->first)
                vertexMap.insert(lb, std::make_pair(polar1, Point(x1, y1))); // angle is not yet in map -> insert
                
            // check wheter angle2 is already in map or not
            lb = vertexMap.lower_bound(polar2);
            if (lb == vertexMap.end() || polar2 != lb->first)
                vertexMap.insert(lb, std::make_pair(polar2, Point(x2, y2))); // angle is not yet in map -> insert
        }
        
        // now transform this map into a sorted list as we can then do very efficient in-place Sutherland-Hodgeman ploygon clipping
        for (std::map<double, Point>::iterator it = vertexMap.begin(); it != vertexMap.end(); ++it)
            outVoronoiRegions[i].push_back(tgt::dvec2(it->second[0], it->second[1]));
        
    }
}

/*** Voreen adaptations end ***/

 
void Delaunay::Triangulate(std::string& triswitches){
    typedef struct triangulateio  TriangStruct;
    typedef struct triangulateio* pTriangStruct;

    in = new TriangStruct;
    pTriangStruct pin = (struct triangulateio *)in;
    
    /*** Voreen adaptations ***/

    vorout = new TriangStruct;
    pTriangStruct pvorout = (struct triangulateio*) vorout;
    
    pin->numberofpoints = (int)PList.size();
    pin->numberofpointattributes = (int)0;
    pin->pointlist = static_cast<double *> ((void *)(&PList[0]));
    pvorout->pointlist = NULL;
    pin->pointattributelist = NULL;         pvorout->pointattributelist = NULL;
    pin->pointmarkerlist = (int *) NULL;    pvorout->pointmarkerlist = (int *) NULL;
    pin->numberofsegments = 0;              pvorout->numberofsegments = 0;
    pin->numberofholes = 0;                 pvorout->numberofholes = 0;
    pin->numberofregions = 0;               pvorout->numberofregions = 0;
    pin->regionlist = (REAL *) NULL;        pvorout->regionlist = (REAL *) NULL;

    pvorout->edgelist = (int*) NULL;
    pvorout->edgemarkerlist = (int*) NULL;
    pvorout->normlist = (REAL*) NULL;

    /*** Voreen adaptations end ***/
    
    delclass = new piyush;
    piyush *pdelclass = (piyush *)delclass;
    triswitches.push_back('\0');
    char *ptris = &triswitches[0];

    pmesh = new piyush::__pmesh;
    pbehavior = new piyush::__pbehavior;

    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
    piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) pbehavior;

    pdelclass->triangleinit(tpmesh);
    pdelclass->parsecommandline(1, &ptris, tpbehavior);

    pdelclass->transfernodes(tpmesh, tpbehavior, pin->pointlist, 
                pin->pointattributelist,
                pin->pointmarkerlist, pin->numberofpoints,
                pin->numberofpointattributes);
    tpmesh->hullsize = pdelclass->delaunay(tpmesh, tpbehavior);

  /* Ensure that no vertex can be mistaken for a triangular bounding */
  /*   box vertex in insertvertex().                                 */
      tpmesh->infvertex1 = (piyush::vertex) NULL;
      tpmesh->infvertex2 = (piyush::vertex) NULL;
      tpmesh->infvertex3 = (piyush::vertex) NULL;

    /* Calculate the number of edges. */
    tpmesh->edges = (3l * tpmesh->triangles.items + tpmesh->hullsize) / 2l;
    
    if (pvorout != (struct triangulateio *) NULL) {
        pvorout->numberofpoints = tpmesh->triangles.items;
        pvorout->numberofpointattributes = tpmesh->nextras;
        pvorout->numberofedges = tpmesh->edges;
    }
    
    pdelclass->numbernodes(tpmesh, tpbehavior);

    /*** Voreen adaptations ***/

    pvorout->edgesOfVoronoiFaces.resize(tpmesh->invertices);
    pdelclass->writevoronoi(tpmesh, tpbehavior, &pvorout->pointlist, &pvorout->pointattributelist,
                            &pvorout->pointmarkerlist, &pvorout->edgelist,
                            &pvorout->edgemarkerlist, &pvorout->normlist, 
                            &pvorout->edgesOfVoronoiFaces[0]);

    /*** Voreen adaptations end ***/                  

    Triangulated = true;
}


Delaunay::~Delaunay(){
        struct triangulateio *pin = (struct triangulateio *)in;
        struct triangulateio *pvorout = (struct triangulateio *)vorout;

    piyush *pdelclass =  (piyush *)delclass;

        piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
        piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) pbehavior;

    pdelclass->triangledeinit(tpmesh, tpbehavior);

    delete tpmesh;
    delete tpbehavior;

    /*** Voreen adaptations ***/
    
    if (pin->pointattributelist) 
        free(pin->pointattributelist);
    if (pin->pointmarkerlist) 
        free(pin->pointmarkerlist);
    if (pin->regionlist) 
        free(pin->regionlist);
    delete pin;
    
    if (pvorout->pointlist) 
        free(pvorout->pointlist);
    if (pvorout->pointattributelist) 
        free(pvorout->pointattributelist);
    if (pvorout->pointmarkerlist)
        free(pvorout->pointmarkerlist);
    if (pvorout->regionlist)
        free(pvorout->regionlist);
    if (pvorout->edgelist) 
        free(pvorout->edgelist);
    if (pvorout->edgemarkerlist) 
        free(pvorout->edgemarkerlist);
    if (pvorout->normlist)
        free(pvorout->normlist);

    delete pvorout;

    /*** Voreen adaptations end ***/ 

    delete pdelclass;    
}

void Delaunay::writeoff(std::string& fname){
    if(!Triangulated) {
        cerr << "FATAL: Write called before triangulation\n";
        exit(1);
    }
        piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
        piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) pbehavior;
    piyush *pdelclass =  (piyush *)delclass;
    char *pfname = new char[fname.size()+1];
    strcpy(pfname , fname.c_str());

    pdelclass->writeoff(tpmesh, tpbehavior, pfname, 0, NULL);
    delete [] pfname;
}

int Delaunay::nedges(){
    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
    return tpmesh->edges;
}

int Delaunay::ntriangles(){
    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
    return tpmesh->triangles.items;    
}

int Delaunay::nvertices(){
    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
        piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) pbehavior;
    int outvertices;

    if (tpbehavior->jettison) {
        outvertices = tpmesh->vertices.items - tpmesh->undeads;
    } else {
        outvertices = tpmesh->vertices.items;
    }

    return outvertices;
}

int Delaunay::hull_size(){
    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
    return  tpmesh->hullsize;    
}

int Delaunay::vertexId(vIterator const &vit){
 piyush::__pmesh     * tpmesh     = (piyush::__pmesh *) vit.MyDelaunay->pmesh;
 return ((int *)vit.vloop)[tpmesh->vertexmarkindex];
}




///////////////////////////////
//
// Vertex Iterator Impl.
//
///////////////////////////////

Delaunay::vIterator::vIterator(Delaunay* adel) {
     //typedef piyush::vertex vertex;
     MyDelaunay = adel;

     piyush::__pmesh     * tpmesh     = (piyush::__pmesh *) adel->pmesh;
     piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) adel->pbehavior;
     piyush *pdelclass =  (piyush *)adel->delclass;

     pdelclass->traversalinit(&( tpmesh->vertices ) );
     vloop = pdelclass->vertextraverse(tpmesh);

     while
        (
         tpbehavior->jettison || 
           (
          ((int *)vloop)[tpmesh->vertexmarkindex+1] == UNDEADVERTEX
        )
    ) 
    vloop = (void *) pdelclass->vertextraverse(tpmesh);
}

Delaunay::vIterator::~vIterator(){
}


Delaunay::vIterator Delaunay::vend(){
    vIterator vit;
    vit.vloop = ((piyush::vertex) NULL);
    return vit;
}



Delaunay::vIterator Delaunay::vIterator::operator++() {
     //typedef piyush::vertex vertex;    

     piyush::__pmesh     * tpmesh     = (piyush::__pmesh *) MyDelaunay->pmesh;
     piyush::__pbehavior * tpbehavior = 
                (piyush::__pbehavior *) MyDelaunay->pbehavior;

     piyush *pdelclass =  (piyush *) MyDelaunay->delclass;

     while
        (
                 tpbehavior->jettison ||
                (
                  ((int *)vloop)[tpmesh->vertexmarkindex+1] == UNDEADVERTEX
                )
        )
        vloop = (void *) pdelclass->vertextraverse(tpmesh);
        vloop = (void *) pdelclass->vertextraverse(tpmesh);

        vIterator vit;
        vit.vloop = vloop;
        vit.MyDelaunay = MyDelaunay;

        return vit;
}

Delaunay::Point & Delaunay::vIterator::operator*() const{
    return *((Point *)vloop);
}


bool operator==(Delaunay::vIterator const &vit1,
            Delaunay::vIterator const &vit2) {
    if (vit1.vloop == vit2.vloop) return true;
    return false;
}


bool operator!=(Delaunay::vIterator const &vit1,
            Delaunay::vIterator const &vit2) {
    if (vit1.vloop != vit2.vloop) return true;
    return false;
}




///////////////////////////////
//
// Face Iterator Impl.
//
///////////////////////////////

Delaunay::fIterator::fIterator(Delaunay* adel) {
     //typedef piyush::vertex vertex;
     typedef piyush::__otriangle trianglelooptype;

     MyDelaunay = adel;

     piyush::__pmesh     * tpmesh     = (piyush::__pmesh *) adel->pmesh;
     piyush *pdelclass =  (piyush *)adel->delclass;

     pdelclass->traversalinit(&( tpmesh->triangles ) );
     // floop = new trianglelooptype;
     trianglelooptype *ploop = (trianglelooptype *)(&floop);
     ploop->tri    = pdelclass->triangletraverse(tpmesh);
     ploop->orient = 0;

}

Delaunay::fIterator::~fIterator(){
}

Delaunay::fIterator Delaunay::fend(){
        fIterator fit;
        //typedef piyush::__otriangle trianglelooptype;
        fit.floop.tri = (double ***) NULL;
        return fit;
}

void Delaunay::fIterator::operator++() {
     // cout << "++ called\n";
     //typedef piyush::vertex   vertex;
     //typedef piyush::triangle triangle;
     typedef piyush::__otriangle trianglelooptype;

     piyush::__pmesh     * tpmesh     = (piyush::__pmesh *) MyDelaunay->pmesh;
     
     trianglelooptype *ploop = (trianglelooptype *)(&floop);
     piyush *pdelclass =  (piyush *) MyDelaunay->delclass;

     ploop->tri = pdelclass->triangletraverse(tpmesh);
     // cout << "tri val = " << ploop->tri << endl;
}

bool operator==(Delaunay::fIterator const &fit1,
                Delaunay::fIterator const &fit2) {

    return (fit1.floop.tri == fit2.floop.tri);

}

bool operator!=(Delaunay::fIterator const &fit1,
                Delaunay::fIterator const &fit2) {
        return !( operator==(fit1,fit2) );
}


/*  A triangle abc has origin (org) a,destination (dest) b, and apex (apex)  */
/*  c.  These vertices occur in counterclockwise order about the triangle.   */
int Delaunay::Org (fIterator const & fit){
     typedef piyush::vertex   vertex;
     //typedef piyush::triangle triangle;
     typedef piyush::__otriangle trianglelooptype;

     piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) 
                        ((fit.MyDelaunay)->pbehavior);

     piyush::__pmesh  * tpmesh  = (piyush::__pmesh *) (fit.MyDelaunay->pmesh);
     trianglelooptype * ploop   = (trianglelooptype *)(&(fit.floop));

     vertex vertexptr = (vertex) ((ploop->tri)[plus1mod3[ploop->orient] + 3]);
     return
    ( ((int *) (vertexptr))[tpmesh->vertexmarkindex] )
    -
    tpbehavior->firstnumber;
        
}


int Delaunay::Dest(fIterator const & fit){
     typedef piyush::vertex   vertex;
     //typedef piyush::triangle triangle;
     typedef piyush::__otriangle trianglelooptype;

     piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *)
                                                ((fit.MyDelaunay)->pbehavior);

     piyush::__pmesh  * tpmesh  = (piyush::__pmesh *) (fit.MyDelaunay->pmesh);
     trianglelooptype * ploop   = (trianglelooptype *)(&(fit.floop));

     vertex vertexptr = (vertex) ((ploop->tri)[minus1mod3[ploop->orient] + 3]);
     return
        ( ((int *) (vertexptr))[tpmesh->vertexmarkindex] )
        -
        tpbehavior->firstnumber;

}

int Delaunay::Apex(fIterator const & fit){
     typedef piyush::vertex   vertex;
     //typedef piyush::triangle triangle;
     typedef piyush::__otriangle trianglelooptype;

     piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *)
                                                ((fit.MyDelaunay)->pbehavior);

     piyush::__pmesh  * tpmesh  = (piyush::__pmesh *) (fit.MyDelaunay->pmesh);
     trianglelooptype * ploop   = (trianglelooptype *)(&(fit.floop));

     vertex vertexptr = (vertex) ((ploop->tri)[ploop->orient + 3]);
     return
        ( ((int *) (vertexptr))[tpmesh->vertexmarkindex] )
        -
        tpbehavior->firstnumber;
}


int Delaunay::Sym(fIterator const & fit, char i){
     typedef piyush::vertex      vertex;
     typedef piyush::triangle    triangle;
     typedef piyush::__otriangle trianglelooptype;
     triangle ptr;                         /* Temporary variable used by sym(). */

     //piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *)
     //                                           ((fit.MyDelaunay)->pbehavior);

     piyush::__pmesh  * tpmesh  = (piyush::__pmesh *) (fit.MyDelaunay->pmesh);
     trianglelooptype * ploop   = (trianglelooptype *)(&(fit.floop));

     char oval = ploop->orient;
     ploop->orient = i;

     trianglelooptype top;
     sym(*ploop, top);
     ploop->orient = oval;

     if(top.tri != tpmesh->dummytri){
         vertex farvertex;
         apex(top, farvertex);
         return ((int *)farvertex)[tpmesh->vertexmarkindex];
     }
     
    return -1;     
}


Delaunay::fIterator Delaunay::Sym(fIterator const & fit){
     fIterator retval;
     retval.MyDelaunay = fit.MyDelaunay;

     //typedef piyush::vertex      vertex;
     typedef piyush::triangle    triangle;
     typedef piyush::__otriangle trianglelooptype;
     triangle ptr;                         /* Temporary variable used by sym(). */

     piyush::__pmesh  * tpmesh  = (piyush::__pmesh *) (fit.MyDelaunay->pmesh);
     trianglelooptype * ploop   = (trianglelooptype *)(&(fit.floop));

     
     trianglelooptype top;
     sym(*ploop, top);
     

     if(top.tri != tpmesh->dummytri){
         retval.floop.tri = top.tri;
         retval.floop.orient = top.orient;

         return retval;
     }
     
     retval.floop.tri = NULL;
     retval.floop.orient = 0;
     return retval;

}

double Delaunay::area(fIterator const & fit){
    Point torg, tdest, tapex;
    torg  = point_at_vertex_id(Org(fit));
    tdest = point_at_vertex_id(Dest(fit));
    tapex = point_at_vertex_id(Apex(fit));
    double dxod(torg[0] - tdest[0]);
    double dyod(torg[1] - tdest[1]);
    double dxda(tdest[0] - tapex[0]);
    double dyda(tdest[1] - tapex[1]);

    double area = 0.5 * (dxod * dyda - dyod * dxda);
    return area;
}

Delaunay::fIterator Delaunay::locate(int vertexid){
    fIterator retval;
    retval.MyDelaunay = this;

    //typedef piyush::vertex      vertex;
    typedef piyush::triangle    triangle;
    typedef piyush::__otriangle trianglelooptype;

    trianglelooptype horiz;               /* Temporary variable for use in locate(). */
    triangle ptr;                         /* Temporary variable used by sym(). */

    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
    piyush::__pbehavior * tpbehavior = (piyush::__pbehavior *) pbehavior;
    piyush  *pdelclass               = (piyush *)              delclass;

    horiz.tri = tpmesh->dummytri;
    horiz.orient = 0;
    symself(horiz);
    double dv[2];
    dv[0] = PList[vertexid][0];
    dv[1] = PList[vertexid][1];

    /* Search for a triangle containing `newvertex'. */
    int intersect = pdelclass->locate(tpmesh, tpbehavior, dv, &horiz);
    if(intersect != piyush::ONVERTEX) { // Not on vertex!
        cout << "Something went wrong in point location\n";
        exit(1);
    }
    
    retval.floop.tri = horiz.tri;
    retval.floop.orient = horiz.orient;

    return retval;

}

inline Delaunay::fIterator Delaunay::Lnext(fIterator const & fit){
/* lnext() finds the next edge (counterclockwise) of a triangle.        */
    typedef piyush::__otriangle trianglelooptype;

    fIterator retval;
    retval.MyDelaunay = this;

    lnext(   (*(trianglelooptype *)(&(fit.floop))), (*(trianglelooptype *)(&(retval.floop))));
    return retval;
}

inline Delaunay::fIterator Delaunay::Lprev(fIterator const & fit){
/*                                                                           */
/*  lprev:  Find the previous edge (clockwise) of a triangle.                */
/*  lprev(abc) -> cab                                                        */
/*                                                                           */
    typedef piyush::__otriangle trianglelooptype;

    fIterator retval;
    retval.MyDelaunay = this;

    lprev(   (*(trianglelooptype *)(&(fit.floop))), (*(trianglelooptype *)(&(retval.floop))));
    return retval;
}


inline Delaunay::fIterator Delaunay::Onext(fIterator const & fit){

/*  onext:  Find the next edge counterclockwise with the same origin.        */
/*  onext(abc) -> ac*                                                        */
    
    typedef piyush::__otriangle trianglelooptype;
    typedef piyush::triangle    triangle;

    triangle ptr;
    fIterator retval;
    retval.MyDelaunay = this;

    //cout << "Onext called:\n " 
    //     << Org(fit) << "\t" << Dest(fit) << "\t" << Apex(fit) << "\n";

    onext(   (*(trianglelooptype *)(&(fit.floop))), (*(trianglelooptype *)(&(retval.floop))));

    // retval could be dummy!
    return retval;
}


inline Delaunay::fIterator Delaunay::Oprev(fIterator const & fit){

/*                                                                           */
/*  oprev:  Find the next edge clockwise with the same origin.               */
/*  oprev(abc) -> a*b                                                        */
/*                                                                           */
    typedef piyush::__otriangle trianglelooptype;
    typedef piyush::triangle    triangle;

    triangle ptr;
    fIterator retval;
    retval.MyDelaunay = this;

    oprev(   (*(trianglelooptype *)(&(fit.floop))), (*(trianglelooptype *)(&(retval.floop))));
    return retval;
}


inline bool Delaunay::isdummy(fIterator const & fit){
    piyush::__pmesh     * tpmesh     = (piyush::__pmesh *)     pmesh;
    typedef piyush::__otriangle trianglelooptype;

    return ( ((trianglelooptype *)(&(fit.floop)))->tri == tpmesh->dummytri );
}

void Delaunay::trianglesAroundVertex(int vertexid, std::vector<int>& ivv){
    fIterator fit = locate(vertexid);
    ivv.clear();

    int start = Dest(fit);
    int linkn = Apex(fit);

    ivv.push_back(vertexid);
    ivv.push_back(start);
    ivv.push_back(linkn);

    fIterator nfit = fit;
    fIterator pnfit = fit; // follows nfit by one triangle

    while( linkn != start ){
        
        nfit = Onext(nfit);
        if (isdummy(nfit)){
            // Do another algorithm
            ivv.clear();

            // use oprev now...
            fit = pnfit;
            nfit = fit;

            start = Apex(fit);
            linkn = Dest(fit);
            
            ivv.push_back(vertexid);
            ivv.push_back(linkn);
            ivv.push_back(start);
            
            while( linkn != start ){
                nfit = Oprev(nfit);
                if(isdummy(nfit))
                    return;
                int a = Org(nfit);
                int b = Dest(nfit);
                int c = Apex(nfit);
                ivv.push_back(a);
                ivv.push_back(b);
                ivv.push_back(c);
                linkn = Dest(nfit);
            }

            return;
        }
        
        pnfit = nfit;

        int a = Org(nfit);
        int b = Dest(nfit);
        int c = Apex(nfit);

        //cout << "Triangle: " << a << "\t"  << b << "\t"  << c << "\n";

        ivv.push_back(a);
        ivv.push_back(b);
        ivv.push_back(c);

        linkn = Apex(nfit);
    }

}

} // namespace tpp ends.
