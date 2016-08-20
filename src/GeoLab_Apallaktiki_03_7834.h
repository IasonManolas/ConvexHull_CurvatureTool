#include "settings.h"
#include "scene.h"
#include "mesh.h"
#include <vector>
#include <string>
#include <MathGeoLib/MathGeoLib.h>
#include <C2DPoint.h>
#include <c:\Users\Iason\Desktop\GeoLab\VVRFramework\include\MathGeoLib\Geometry\Triangle.h>
using namespace vvr;

class qhVertex
{
public:
	std::shared_ptr<class qhHalfEdge> Edge;// one of the half-edges emantating from the vertex
	vvr::Vec3d Position;
	vector<std::shared_ptr<class qhFace>> ConflictList;
	qhVertex(){}
	qhVertex(vvr::Vec3d Position3D){
		Position = Position3D;
	}
	qhVertex(std::shared_ptr<class qhHalfEdge> emanating, vvr::Vec3d Position3D){
		Edge = emanating;
		Position = Position3D;
	}

};

class qhHalfEdge
{
public:
	std::shared_ptr<class qhVertex> Begin;// vertex at the BEGINNING of the half-edge
	std::shared_ptr<class qhHalfEdge> Prev;
	std::shared_ptr<qhHalfEdge> Next;// next half-edge around the face
	std::shared_ptr<qhHalfEdge> Twin;// oppositely oriented adjacent half-edge 
	std::shared_ptr<qhFace> Face;// face the half-edge borders
	qhHalfEdge(){}
	qhHalfEdge(std::shared_ptr<qhVertex> begin)
	{
		Begin = begin;
	}
	qhHalfEdge(std::shared_ptr<class qhVertex> begin, std::shared_ptr<class qhHalfEdge> twin, std::shared_ptr<class qhFace> face, std::shared_ptr<class qhHalfEdge> next, std::shared_ptr<class qhHalfEdge> previous)
	{
		Begin = begin;
		Prev = previous;
		Next = next;
		Twin = twin;
		Face = face;
	}
	std::shared_ptr<qhHalfEdge> findPrevious()
	{
		return this->Next->Next;
	}
	void changeRotation() // does not change twins. To change twins : edge->twin->face->changeRotation() for all edges of face
	{
		std::shared_ptr<qhHalfEdge> NewNext, NewPrev,NewTwin;
		std::shared_ptr<class qhVertex> NewBegin;

		NewNext = this->Prev;
		NewBegin = this->Next->Begin;
		NewPrev = this->Next;
		
		//face stays the same

		this->Begin = NewBegin;
		this->Next = NewNext;
		this->Prev = NewPrev;

	}
};
class qhFace
{
public:
	std::shared_ptr<class qhHalfEdge> Edge;
	vector<std::shared_ptr<class qhVertex>> ConflictList;
	qhFace(){}
	qhFace(std::shared_ptr<class qhHalfEdge> edge)
	{
		Edge = edge;
	}
	vvr::Vec3d CMofFace()
	{

		double cmx = (this->Edge->Begin->Position.x + this->Edge->Next->Begin->Position.x + this->Edge->findPrevious()->Begin->Position.x) / 3,
			cmy = (this->Edge->Begin->Position.y + this->Edge->Next->Begin->Position.y + this->Edge->findPrevious()->Begin->Position.y) / 3,
			cmz = (this->Edge->Begin->Position.z + this->Edge->Next->Begin->Position.z + this->Edge->findPrevious()->Begin->Position.z) / 3;
		return (cmx, cmy, cmz);

	}
	void triChangeRotation()
	{
		std::shared_ptr<qhHalfEdge> e1, e2, e3;
		e1 = this->Edge;
		std::shared_ptr<class qhVertex> e3Begin = this->Edge->Begin;
		e2 = e1->Next;
		e3 = e1->Prev;
		e1->changeRotation();
		e2->changeRotation(); 
		e3->changeRotation();
		e3->Begin = e3Begin;
		e1->Twin->Begin;
	}
};
struct point2D
{
	int indexInPointsToProject;
	C2DPoint coordinates2D;
};
struct vertex{
	int indexInVerts;
	double verKmCurv;
};
class Simple3DScene : public vvr::Scene 
{
public:
    Simple3DScene();
    const char* getName() const { return "Simple 3D Scene";}
    void keyEvent(unsigned char key, bool up, int modif) override;
    void mousePressed(int x, int y, int modif) override;
    void mouseMoved(int x, int y, int modif) override;

protected:
    void draw() override;
    void resize() override;

private:
	//Erwthma 1
	void Curvature2D();
    void drawPolygon2D();
	void Curvature3D();
	void drawPolygon3D();
    void savePolygonToFile();
    void loadPolygonFromFile(string filename);
    void pixelCoordsToSceneCoords(float &x, float &y);
	//Erwthma 2
	void CurvOf3DModel();
	void oneRingNeighbours();
	void neighbourTriangles();
	void normalsOfvertices();
	void neighbourVertices();
	void tiAndKnOfTi();
	void abcCoefCacl();
	void Kg();
	void Km();
	void curvOfTris();
	void DrawKm();
	void DrawKg();
	//Erwthma 3
	void ConvexHull3D();
	void CHInitialization();
	void CHExpand(int);
	void CHConflicts(int);
	void CHInitialConflicts();
	void PiHorizon(int, vector<vector<int>> &);
	void DrawHorizon(vector<vector<int>>);
	void CreateNewFaces(int, vector<vector<int>>);
	void getHorizonEdges(vector<vector<int>> adjEdges, int triIndex, vector<vector<int>> &horizonEdges);
	math::Triangle CHToMathTr(int);
	void newCM(int);
	bool IsCCW(int Tri[3]);
	void DrawCH();
	void DeleteVisibleFaces(int);
	void DrawConflictsOf(int);
	void curvOfConvexHull();
	void KmCH();
	void DrawKmCH();
	//Erwthma 4
	void DrawRetriKmCH();
	void Retriangulation_Of_CH();
	void CreateFourTriangles(Vec3d a, Vec3d b, Vec3d c, Vec3d mid1, Vec3d mid2, Vec3d mid3, vector<Vec3d> &Tr1, vector<Vec3d>& Tr2, vector<Vec3d> &Tr3, vector<Vec3d> &Tr4);
	math::Triangle triCHtoMathTr(int);
	void Curv_Of_Retri_CH();
	void ReRetriangulate();
	//Erwthma 5
	void Topika_Megista();
	bool IsLocalExtrema(vertex);
	void DrawExtremaPoints(int);
	//Erwthma 6
	void k_Daktylios(int, double);
	int megethos_Daktyliou();
	double katwfli_Kampylothtas();
	void diwrthwsh_ena_daktyliou(double katwfli);
	void diwrthwsh_dyo_daktyliou(double katwfli);
	void diwrthwsh_tria_daktyliou(double katwfli);
	void curvOfPoint(int);
	void normalsOfvertices(int index);
	void tiAndKnOfTi(int index);
	void abcCoefCacl(int index);
	void Km(int index);
	void Kg(int index);
	void twoRingNeighbours();
	void threeRingNeighbours();
	void twoRingExtremePoints();
	void threeRingExtremePoints();
	bool isTwoRingLocalExtr(int index);
	bool isThreeRingLocalExtr(int index);
	void drawTwoRingLocalExtrPoints();
	void drawThreeRingLocalExtrPoints();
private:
    bool                    b_show_pts;
    vector<C2DPoint>        m_pts;
	vector<C3DPoint>        m_pts3D;
    vvr::Mesh               m_model;
	vvr::Mesh				m_convexHull;
    vvr::Settings           m_settings;
    vvr::Colour             m_obj_col;
    float                   m_sphere_rad;
    int                     m_style_flag;
	vector<vector<int>>		neigTri; //exei megethos oses oi koryfes . Kathe grammh i exei ena vector pou deixnei se poia trigwna yparxei h i koryfh.
	vector<vvr::Vec3d>		verNorm; //exei megethos oses oi koryfes . Kathe grammh i deixnei to normal dianysma sthn vertex[i] .
	vector<vector<int>>		verNeig; //exei gia kathe koryfh i ena vector me ta index twn geitonwn tou i sto vector mVertices. 
	vector<vector<vvr::Vec3d>> verTi; // Vriskei gia kathe koryfh i ta dianysmata ti ths koryfhs i me olous tous geitones ths ta vazei se vector kai ayton ton vazei sth thesh verTi[i].
	vector<vector<double>>	knofti; // Exei megethos oso to verts. Gia kathe koryfh i vriskei thn kampylothta kn(ti) gia kathe geitona ths i tis vazei se vector kai ayton sth thesh knofti
	vector<vector<double>>  verABC; //Exei gia kathe koryfh i ena vector triwn stoixeiwn tou opoiou to prwto stoixeio einai to kn(tid) to deutero to b kai to trito to c.
	vector<double>			verKg; //Vector ston opoio to i stoixeio einai h Gaussian kampylothta ths i koryfhs.
	vector<double>			verKm;	//Vector ston opoio to i stoixeio einai h Messh kampylothta ths i koryfhs.
	vector<double>			triKg;
	vector<double>			triKm;
	
	vector<vector<int>>		CH;
	vector<int>				P;//Vector of points tha havent yet been prosecced.
	vector<vector<int>>		Fconf; // Fconf(pi) ,witch faces does pi see ,pi: point that hasnt been added to CH .
	vector<vector<int>>		Pconf; //Pconf(fi) ,  fi: face , witch points can see fi.
	vector<int>				CHVerts;
	Vec3d					CHCM; // Convex Hulls Centre of Mass.

	vector<vector<Vec3d>>	vertsOfTriCH;
	vector<double>			triKmCH;
	vector<vector<Vec3d>>	triangulatedCH;
	vector<double>			reTriKmCH;
	vector<vertex>			vertsIndexAndCurv; // Periexei oles tis koryfes. Exei to index ths koryfhs kai thn messh kampylothta ths.
	vector<vertex>			extremaPoints;
	vector<Vec3d>			problematicPoints;
	vector<vector<int>>		twoRingNeig;
	vector<vector<int>>		threeRingNeig;
	vector<int>			twoRingExtrPoints;
	vector<int>			threeRingExtrPoints;
};
class 
Vec3d crossProduct(Vec3d a, Vec3d b)
{
	Vec3d CP;
	CP.x = (a.y*b.z - a.z*b.y);
	CP.y = (a.z*b.x - a.x*b.z);
	CP.z=(a.x*b.y - a.y*b.x);
	return CP;
}
double dotProduct(Vec3d a, Vec3d b)
{
	double x, y, z, DP;
	x = a.x*b.x;
	y = a.y*b.y;
	z = a.z*b.z;
	DP = x + y + z;
	return DP;
}
Vec3d subVec(Vec3d a, Vec3d b)
{
	Vec3d vec;
	vec.x = a.x - b.x;
	vec.y = a.y - b.y;
	vec.z = a.z - b.z;
	return vec;
}
Vec3d addVec(Vec3d a, Vec3d b)
{
	Vec3d vec;
	vec.x = a.x + b.x;
	vec.y = a.y + b.y;
	vec.z = a.z + b.z;
	return vec;
}
Vec3d scaleVec(Vec3d a, double c)
{
	Vec3d vec;
	vec.x = a.x*c;
	vec.y = a.y*c;
	vec.z = a.z*c;
	return vec;
}
double lengthVec(Vec3d a)
{
	return sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
}
Vec3d norVec(Vec3d a)
{
	return scaleVec(a, 1 / lengthVec(a));
}

vec vec3dToVec(Vec3d a)
{
	vec b;
	b.Set(a.x, a.y, a.z);
	return b;
}
bool Simple3DScene::IsCCW(int Tri[3])
{
	vector<Vec3d> &verts = m_model.getVertices();
	Vec3d p1 = verts[Tri[0]],p2=verts[Tri[1]],p3=verts[Tri[2]];
	vec v1, v2, v3;
	v1 = vec3dToVec(p1);
	v2 = vec3dToVec(p2);
	v3 = vec3dToVec(p3);
	Plane p(v1, v2, v3);
	if (p.IsOnPositiveSide(vec3dToVec(CHCM))) return false;
	return true;
}

math::Triangle Simple3DScene::CHToMathTr(int indexInCH)
{
	vector<Vec3d> &verts = m_model.getVertices();
	vector<int> vec = CH[indexInCH];
	math::Triangle Tri;
	Tri.a.Set(verts[vec[0]].x, verts[vec[0]].y, verts[vec[0]].z);
	Tri.b.Set(verts[vec[1]].x, verts[vec[1]].y, verts[vec[1]].z);
	Tri.c.Set(verts[vec[2]].x, verts[vec[2]].y, verts[vec[2]].z);
	return Tri;
}
