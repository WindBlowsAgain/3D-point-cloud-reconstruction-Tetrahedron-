#ifndef C3D_T_H
#define C3D_T_H
#include"Data_type.h"
#include <vector>
using namespace std;


struct super_tetrahedron_vertices
{
	C3DPoint v1;
	C3DPoint v2;
	C3DPoint v3;
	C3DPoint v4;

};

struct distance_vertices
{
	double d;
	int index;
	bool friend operator<(const distance_vertices &C1, const distance_vertices &C2)
	{
		if (C1.d<C2.d) return true;
		if (C1.d == C2.d && C1.index<C2.index)return true;
		return false;
	}

};
struct angle_vertices
{
	double angle;
	int index;
	bool friend operator<(const angle_vertices &C1, const angle_vertices &C2)
	{
		if (C1.angle<C2.angle) return true;
		if (C1.angle == C2.angle && C1.index<C2.index)return true;
		return false;
	}

};


struct line_triangles
{
	int a, b;//线段两个个点的索引
	vector<int> v_index; //线段所在三角面的索引
	bool friend operator<(const line_triangles  &t1, const line_triangles  &t2)
	{
		if (t1.a<t2.a) return true;
		if (t1.a == t2.a && t1.b<t2.b)return true;
		if (t1.a == t2.a && t1.b == t2.b&&t1.v_index[0]<t2.v_index[0])return true;
		return false;
	}
};
struct triangle_one_vertice
{
	int triangle, vertice;//线段两个个点的索引
	double angle;
	bool friend operator<(const triangle_one_vertice  &t1, const triangle_one_vertice  &t2)
	{
		if ((-t1.angle)<(-t2.angle)) return true;
		if ((-t1.angle) == (-t2.angle) && t1.vertice<t2.vertice)return true;
		return false;
	}

};

class C3Delaunay_tetrahedron
{
public:
	vector<C3DPoint> vertices;
	C3Delaunay_tetrahedron(vector<C3DPoint>);
	super_tetrahedron_vertices super_tetrahedron();
	tetrahedron tetrahedron_circumsphere(int, int, int, int);
	vector<C3DPoint> mesh_generation();
	vector<C3DPoint>  mesh_triangles(int, vector<int>);//以核心点为中心生成局部小三角形



};

#endif