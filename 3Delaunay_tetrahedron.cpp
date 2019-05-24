#include "stdafx.h"
#include "3Delaunay_tetrahedron.h"
#include <vector>
using namespace std;
C3Delaunay_tetrahedron::C3Delaunay_tetrahedron(vector<C3DPoint> temp_vertices)
{
	vertices.insert(vertices.begin(), temp_vertices.begin(), temp_vertices.end());
}
super_tetrahedron_vertices C3Delaunay_tetrahedron::super_tetrahedron()//���������Ľ���
{
	double 	xmin, xmax, ymin, ymax, dx, dy, dmax,zmin,zmax,dz;
	int v_size = vertices.size();
	xmin = vertices[0].x;
	xmax = vertices[v_size - 1].x;
	ymin = vertices[0].y;
	ymax = vertices[0].y;
	zmin = vertices[0].z;
	zmax = vertices[0].z;
	int i = 0;
	for (i = 0; i < v_size; i++)
	{
		if (vertices[i].y >= ymax) ymax = vertices[i].y;
		if (vertices[i].y <= ymin) ymin = vertices[i].y;
		if (vertices[i].z >= zmax) zmax = vertices[i].z;
		if (vertices[i].z <= zmin) zmin = vertices[i].z;
	}
	
	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;
	xmin = xmin - 3*dx;//ȷ���㼯��û�е����ڳ��������� ��
	ymin = ymin - 3*dy ;
	zmin = zmin - 3*dz ;


	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;
	super_tetrahedron_vertices temp;
	//P0(Xmin��Ymin��Zmin)
	temp.v1.x = xmin;
	temp.v1.y = ymin; 
	temp.v1.z = zmin;
	//P1(Xmin��Ymin��Zmin + 4 * dz)
	temp.v2.x = xmin;
	temp.v2.y = ymin;
	temp.v2.z = zmin+4*dz;
	//P2(Xmin��Ymin+ 4* dy��Zmin)
	temp.v3.x = xmin;
	temp.v3.y = ymin+4*dy;
	temp.v3.z = zmin;
	//P3(Xmin+ 4*dx��Ymin��Zmin)
	temp.v4.x = xmin+4*dx;
	temp.v4.y = ymin;
	temp.v4.z = zmin;
	return temp;
}
tetrahedron C3Delaunay_tetrahedron::tetrahedron_circumsphere(int i, int j, int k, int w)////���������
{

	double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
	x1 = vertices[i].x;
	y1 = vertices[i].y;
	z1 = vertices[i].z;
	x2 = vertices[j].x;
	y2 = vertices[j].y;
	z2 = vertices[j].z;
	x3 = vertices[k].x;
	y3 = vertices[k].y;
	z3 = vertices[k].z;
	x4 = vertices[w].x;
	y4 = vertices[w].y;
	z4 = vertices[w].z;
	double a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3, d, d1, d2, d3,x,y,z,r2;
	a11 = 2 * (x2 - x1); a12 = 2 * (y2 - y1); a13 = 2 * (z2 - z1);
	a21 = 2 * (x3 - x2); a22 = 2 * (y3 - y2); a23 = 2 * (z3 - z2);
	a31 = 2 * (x4 - x3); a32 = 2 * (y4 - y3); a33 = 2 * (z4 - z3);
	b1 = x2*x2 - x1*x1 + y2*y2 - y1*y1 + z2*z2 - z1*z1;
	b2 = x3*x3 - x2*x2 + y3*y3 - y2*y2 + z3*z3 - z2*z2;
	b3 = x4*x4 - x3*x3 + y4*y4 - y3*y3 + z4*z4 - z3*z3;
	d = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
	d1 = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - b1*a23*a32 - a12*b2*a33 - a13*a22*b3;
	d2 = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - a11*a23*b3 - b1*a21*a33 - a13*b2*a31;
	d3 = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - a11*b2*a32 - a12*a21*b3 - b1*a22*a31;
	x = d1 / d;
	y = d2 / d;
	z = d3 / d;
	r2 = (x1 - x)*(x1 - x) + (y1 - y)*(y1 - y) + (z1 - z)*(z1 - z);

	
	tetrahedron t_temp;
	t_temp.i = i;
	t_temp.j = j;
	t_temp.k = k;
	t_temp.w = w;
	t_temp.xc = x;
	t_temp.yc = y;
	t_temp.zc = z;
	t_temp.r = sqrt(r2);

	return t_temp;
}
vector<C3DPoint>  C3Delaunay_tetrahedron::mesh_triangles(int i, vector<int> v)//�Ժ��ĵ�Ϊ�������ɾֲ�С������
{
	//������������ĵ㣬�������������������������ɿռ���
	//�Ȱ�����������ĵ����������γ������ߣ��Ե�һ������ĵ������Ϊ��׼�ߣ��û�׼�߲�������������ٵ�˿ռ�����ȷ���������ڵĲ��棬�õ��ȷ���нǣ��ۺ�����������ȷ�����������ڵ�λ�ã������������
	/*************************************************///������������
	//�Ե�һ��Ϊ��׼��Ѱ����������������������ĵ㣬���ɻ�׼�ߣ�Ȼ������ȷ�������ĵ����׼�ߵļнǣ���ȷ���ĸ������ڵ�λ�ã������������
	vector<C3DPoint> temp_triangles;
	
		C3DPoint first_vertice(vertices[v[0]]);
		set<distance_vertices> s_distance;
		for (int j = 1; j < v.size(); j++)
		{
			distance_vertices temp;
			temp.d = sqrt((first_vertice.x - vertices[v[j]].x)*(first_vertice.x - vertices[v[j]].x) + (first_vertice.y - vertices[v[j]].y)*(first_vertice.y - vertices[v[j]].y) + (first_vertice.z - vertices[v[j]].z)*(first_vertice.z - vertices[v[j]].z));
			temp.index = v[j];
			s_distance.insert(temp);
		}
		vector<distance_vertices> v_distance;
		v_distance.resize(v.size()-1);
		copy(s_distance.begin(), s_distance.end(), v_distance.begin());
		C3DPoint second_vertice(vertices[v_distance[0].index]);
		//Ѱ����������������������ĵ㣬���ɻ�׼�ߣ�Ȼ������ȷ�������ĵ����׼�ߵļнǣ���ȷ���ĸ������ڵ�λ�ã������������
		set <angle_vertices> s_angle_vertices;//�����׼�ߵļнǽ�������
		vector <angle_vertices> v_angle_vertices;
		vector<int> order_vertices;
		order_vertices.resize(v.size());
		order_vertices[0] = v[0];
		order_vertices[1] = v_distance[0].index;
		for (int j = 1; j < v_distance.size(); j++)
		{
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;
			x1 = first_vertice.x;
			y1 = first_vertice.y;
			z1 = first_vertice.z;
			x2 = second_vertice.x;
			y2 = second_vertice.y;
			z2 = second_vertice.z;
			x3 = vertices[v_distance[j].index].x;
			y3 = vertices[v_distance[j].index].y;
			z3 = vertices[v_distance[j].index].z;
			double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
			double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3) + (z1 - z3)*(z1 - z3));
			double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3) + (z2 - z3)*(z2 - z3));
			//double A = (acos((b*b + c*c - a*a) / (2 * b*c)) / 3.1415926) * 180;
			//double B = (acos((a*a + c*c - b*b) / (2 * a*c)) / 3.1415926) * 180;
			double temp_angle = (acos((a*a + b*b - c*c) / (2 * a*b)) / 3.1415926) * 180;
			
			angle_vertices temp_angle_vertices;
			temp_angle_vertices.angle = temp_angle;
			temp_angle_vertices.index = v_distance[j].index;
			s_angle_vertices.insert(temp_angle_vertices);
		}
		v_angle_vertices.resize(v.size()-2);
		copy(s_angle_vertices.begin(), s_angle_vertices.end(), v_angle_vertices.begin());
		for (int j = 0; j < v_angle_vertices.size(); j++)
		{
			order_vertices[j+2] = v_angle_vertices[j].index;
		}
		for (int k = 0; k < v.size() - 1; k++)
		{
			C3DPoint temp;
			temp.x = i;
			temp.y = order_vertices[k];
			temp.z = order_vertices[k+1];
			temp_triangles.push_back(temp);
		}
		C3DPoint t_emp;
		t_emp.x = i;
		t_emp.y = order_vertices[0];
		t_emp.z = order_vertices[v.size() - 1];
		temp_triangles.push_back(t_emp);

	
	return temp_triangles;
}
vector<C3DPoint> C3Delaunay_tetrahedron::mesh_generation()
{
	//��������������붥������
	int n = vertices.size();
	super_tetrahedron_vertices temp;
	temp = super_tetrahedron();
	C3DPoint p;
	C3DPoint q;
	C3DPoint w;
	C3DPoint f;
	p.x = temp.v1.x;
	p.y = temp.v1.y;
	p.z = temp.v1.z;
	q.x = temp.v2.x;
	q.y = temp.v2.y;
	q.z = temp.v2.z;
	w.x = temp.v3.x;
	w.y = temp.v3.y;
	w.z = temp.v3.z;
	f.x = temp.v4.x;
	f.y = temp.v4.y;
	f.z = temp.v4.z;
	vertices.push_back(p);
	vertices.push_back(q);
	vertices.push_back(w);
	vertices.push_back(f);


	vector<tetrahedron> temp_tetrahedron;//�������ĸ��������������������ģ��뾶
	vector<tetrahedron> final_tetrahedron;
	vector<int> triangles;
	tetrahedron temp_t;
	temp_t = tetrahedron_circumsphere(n, n + 1, n + 2,n+3);
	temp_tetrahedron.push_back(temp_t);

	for (int i = 0; i < n; i++)//�����ؽ�ÿ��������ӵ�������
	{

		for (int j = 0; j < temp_tetrahedron.size(); j++)//����temp triangles�е�ÿһ��������
		{
			double dx, dy,dz;
			dx = vertices[i].x - temp_tetrahedron[j].xc;
			//��������������������������ұߣ���ô���������Ͳ����ٱ�����ˡ��Ӵ��б���ɾ������������ӵ��ر��б��У���������
			if (dx > 0.0 && dx > temp_tetrahedron[j].r)
			{
				final_tetrahedron.push_back(temp_tetrahedron[j]);
				temp_tetrahedron.erase(temp_tetrahedron.begin() + j);
				j--;
				continue;
			}
			//����������������⣬�������������
			dy = vertices[i].y - temp_tetrahedron[j].yc;
			dz = vertices[i].z - temp_tetrahedron[j].zc;
			if (sqrt(dx * dx + dy * dy+dz*dz)>temp_tetrahedron[j].r)
				continue;
			//����������������ڣ����Ƴ������岢��������������ӵ��������б���
			triangles.push_back(temp_tetrahedron[j].i);
			triangles.push_back(temp_tetrahedron[j].j);
			triangles.push_back(temp_tetrahedron[j].k);
			triangles.push_back(temp_tetrahedron[j].i);
			triangles.push_back(temp_tetrahedron[j].j);
			triangles.push_back(temp_tetrahedron[j].w);
			triangles.push_back(temp_tetrahedron[j].i);
			triangles.push_back(temp_tetrahedron[j].k);
			triangles.push_back(temp_tetrahedron[j].w);
			triangles.push_back(temp_tetrahedron[j].j);
			triangles.push_back(temp_tetrahedron[j].k);
			triangles.push_back(temp_tetrahedron[j].w);
			temp_tetrahedron.erase(temp_tetrahedron.begin() + j);
			j--;//��Ϊɾ��һ��Ԫ�غ�ʣ�µĻ�ǰ��һλ
		}

		//ɾ���ظ������Σ�ע�Ⲣ��ֻ�ǽ��ظ��Ĳ���ȥ��, ���ǽ������ظ���trianglesȫ��ȥ��, ������ȥ���ظ��Ĳ���, ԭ����triangles ҲҪȥ��).
		for (int i = 0; i < triangles.size(); i += 3)
		{
			int a = triangles[i];
			int b = triangles[i + 1];
			int c = triangles[i + 2];
			for (int j = i + 3; j< triangles.size(); j += 3)
			{
				int m = triangles[j];
				int n = triangles[j + 1];
				int v = triangles[j + 2];

				if ((a == m && b == n && c == v) || (a == m && c == n && b == v) || (b == m && a == n &&c == v) || (b == m && c == n &&a == v) || (c == m && a == n &&b == v) || (c == m && b == n &&a == v))
				{
					triangles.erase(triangles.begin() + i, triangles.begin() + (i + 3));
					triangles.erase(triangles.begin() + (j - 3), triangles.begin() + j);
					i -= 3;//��Ϊɾ������Ԫ�غ�ʣ�µĻ�ǰ����λ
					break;
				}
			}
		}
		//Ϊÿ�����������һ���µ�������
		for (int j = 0; j < triangles.size(); j += 3)
		{
			int a = triangles[j];
			int b = triangles[j + 1];
			int c = triangles[j + 2];
			temp_tetrahedron.push_back(tetrahedron_circumsphere(a, b,c, i));
		}
		vector<int>().swap(triangles);
	}
	//ɾ���볬�����干��һ������������壬����һ����ʾ������������б�
	for (int i = 0; i < temp_tetrahedron.size(); i++)
	{
		final_tetrahedron.push_back(temp_tetrahedron[i]);
	}

	vector<tetrahedron>v_tetrahedron;//�����������ĸ�����±�
	for (int j = 0; j < final_tetrahedron.size(); j++)
	if (final_tetrahedron[j].i < n && final_tetrahedron[j].j < n &&final_tetrahedron[j].k < n&&final_tetrahedron[j].w < n)
	{
		v_tetrahedron.push_back(final_tetrahedron[j]);
	}
	/***********************************************************************************/

	/***********************************************************************************/


	vector<C3DPoint>v_triangles;//���������������������Ƭ
	set<C3DPoint>s_triangles;
	for (int n = 0; n < v_tetrahedron.size(); n++)
	{
		set<int> s_temp;
		s_temp.insert(v_tetrahedron[n].i);
		s_temp.insert(v_tetrahedron[n].j);
		s_temp.insert(v_tetrahedron[n].k);
		s_temp.insert(v_tetrahedron[n].w);
		vector<int> v_temp;
		v_temp.resize(s_temp.size());
		copy(s_temp.begin(), s_temp.end(), v_temp.begin());
		C3DPoint p;
		C3DPoint q;
		C3DPoint w;
		C3DPoint f;
		p.x = v_temp[0];
		p.y = v_temp[1];
		p.z = v_temp[2];
		q.x = v_temp[0];
		q.y = v_temp[1];
		q.z = v_temp[3];
		w.x = v_temp[0];
		w.y = v_temp[2];
		w.z = v_temp[3];
		f.x = v_temp[1];
		f.y = v_temp[2];
		f.z = v_temp[3];
		s_triangles.insert(p);
		s_triangles.insert(q);
		s_triangles.insert(w);
		s_triangles.insert(f);
		/*v_triangles.push_back(p);
		v_triangles.push_back(q);
		v_triangles.push_back(w);
		v_triangles.push_back(f);*/
	}

	
	v_triangles.resize(s_triangles.size());
	copy(s_triangles.begin(), s_triangles.end(), v_triangles.begin());

	/***********************************************************************************/
	//���ж�������Ƭ�����Բ���ų��ϴ��������

	//һ�ι�������Ƭ
	vector<double> v_r;//�����Ӧ�±�����Ƭ�İ뾶
	double sum_triangles_r = 0;
	double average_triangles_r = 0;
	for (int i = 0; i < v_triangles.size(); i++)
	{
		double x1, y1, z1, x2, y2, z2, x3, y3, z3;
		x1 = vertices[v_triangles[i].x].x;
		y1 = vertices[v_triangles[i].x].y;
		z1 = vertices[v_triangles[i].x].z;
		x2 = vertices[v_triangles[i].y].x;
		y2 = vertices[v_triangles[i].y].y;
		z2 = vertices[v_triangles[i].y].z;
		x3 = vertices[v_triangles[i].z].x;
		y3 = vertices[v_triangles[i].z].y;
		z3 = vertices[v_triangles[i].z].z;
		double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
		double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3) + (z1 - z3)*(z1 - z3));
		double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3) + (z2 - z3)*(z2 - z3));
		double p = (a + b + c) / 2;
		double S = sqrt(p*(p - a)*(p - b)*(p - c));
		double radius = a*b*c / (4 * S);
		sum_triangles_r += radius;
		v_r.push_back(radius);
	}

	average_triangles_r = sum_triangles_r / (v_r.size());
	for (int i = 0; i < v_triangles.size(); i++)
	{
		if (v_r[i]>average_triangles_r)
		{
			v_triangles.erase(v_triangles.begin() + i);
			v_r.erase(v_r.begin() + i);
			i--;
		}
	}

	/***********************************************************************************/
	//���ι�������Ƭ
	vector<double> v_r_2;//�����Ӧ�±�����Ƭ�İ뾶
	double sum_triangles_r_2 = 0;
	double average_triangles_r_2 = 0;
	for (int i = 0; i < v_triangles.size(); i++)
	{
		double x1, y1, z1, x2, y2, z2, x3, y3, z3;
		x1 = vertices[v_triangles[i].x].x;
		y1 = vertices[v_triangles[i].x].y;
		z1 = vertices[v_triangles[i].x].z;
		x2 = vertices[v_triangles[i].y].x;
		y2 = vertices[v_triangles[i].y].y;
		z2 = vertices[v_triangles[i].y].z;
		x3 = vertices[v_triangles[i].z].x;
		y3 = vertices[v_triangles[i].z].y;
		z3 = vertices[v_triangles[i].z].z;
		double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
		double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3) + (z1 - z3)*(z1 - z3));
		double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3) + (z2 - z3)*(z2 - z3));
		double p = (a + b + c) / 2;
		double S = sqrt(p*(p - a)*(p - b)*(p - c));
		double radius = a*b*c / (4 * S);
		sum_triangles_r_2 += radius;
		v_r_2.push_back(radius);
	}

	average_triangles_r_2 = sum_triangles_r_2 / (v_r_2.size());
	for (int i = 0; i < v_triangles.size(); i++)
	{
		if (v_r_2[i]>average_triangles_r_2*1.3)
		{
			v_triangles.erase(v_triangles.begin() + i);
			v_r_2.erase(v_r_2.begin() + i);
			i--;
		}
	}
	
	/***********************************************************************************/
//    ��һ�����Ӧ���߶��γ�һ�����ϣ�ȡ��Щ�߶��о�����̵������ֱ��γ�������
	vector<vector<int>>   v_vertices_triangles;//�����Ӧ�±�����ڵ������Σ����������ɵģ�
	vector<set<int>>   v_vertices_linked_vertices;//�������Ӧ�±���Ѿ����ߵĵ���������������������Σ�
	vector<C3DPoint>  v_output_triangles;//�����������ɵ�������
	v_vertices_triangles.resize(vertices.size()-4);
	v_vertices_linked_vertices.resize(vertices.size() - 4);
	
	for (int i = 0; i < v_triangles.size(); i++)
	{
		v_vertices_triangles[v_triangles[i].x].push_back(i);
		v_vertices_triangles[v_triangles[i].y].push_back(i);
		v_vertices_triangles[v_triangles[i].z].push_back(i);

	}
	int t = 0;
	for (int i = 0; i < v_vertices_triangles.size(); i++)
	{
		vector<int>  v_vertices_mindistance_vertices;//�������Ӧ�±������������������
		set<int> s_vertices;
		vector<int> v_vertices;
		vector<int>linked_vertices;
		if (v_vertices_linked_vertices[i].size()>6) continue;
		if (v_vertices_triangles[i].size() == 0)
		{
			t++;
			continue;
		}//���ĳ�㲻����κ������Σ���Ϊ������

		for (int j = 0; j < v_vertices_triangles[i].size(); j++)
		{
			if (v_vertices_linked_vertices[v_triangles[(v_vertices_triangles[i])[j]].x].size()<9)
			s_vertices.insert(v_triangles[(v_vertices_triangles[i])[j]].x);
			if (v_vertices_linked_vertices[v_triangles[(v_vertices_triangles[i])[j]].y].size()<9)
			s_vertices.insert(v_triangles[(v_vertices_triangles[i])[j]].y);
			if (v_vertices_linked_vertices[v_triangles[(v_vertices_triangles[i])[j]].z].size()<9)
			s_vertices.insert(v_triangles[(v_vertices_triangles[i])[j]].z);
			
			/*s_vertices.insert(v_triangles[(v_vertices_triangles[i])[j]].x);
			
			s_vertices.insert(v_triangles[(v_vertices_triangles[i])[j]].y);
			
			s_vertices.insert(v_triangles[(v_vertices_triangles[i])[j]].z);*/
		}
		s_vertices.erase(i);
		v_vertices.resize(s_vertices.size());
		copy(s_vertices.begin(), s_vertices.end(), v_vertices.begin());


	

		if (s_vertices.size() <= 6)
		{
			v_vertices_mindistance_vertices.resize(s_vertices.size());
			copy(s_vertices.begin(), s_vertices.end(), v_vertices_mindistance_vertices.begin());
		}
		else
		{

			set<distance_vertices> s_distance;
			
			for (int j = 0; j < v_vertices.size(); j++)
			{
				vector<int> v_linked_vertices_temp;//�����жϸõ��Ƿ��Ѿ���i������
				v_linked_vertices_temp.resize(v_vertices_linked_vertices[i].size());
				copy(v_vertices_linked_vertices[i].begin(), v_vertices_linked_vertices[i].end(), v_linked_vertices_temp.begin());
				
				bool temp = FALSE;
				for (int k = 0; k < v_linked_vertices_temp.size(); k++)
				{
					if (v_vertices[j] == v_linked_vertices_temp[k])
						temp = TRUE;

				}
				if (temp == TRUE)
				{
					linked_vertices.push_back(v_vertices[j]);
				}
				else
				{
					distance_vertices temp;
					temp.d = sqrt((vertices[i].x - vertices[v_vertices[j]].x)*(vertices[i].x - vertices[v_vertices[j]].x) + (vertices[i].y - vertices[v_vertices[j]].y)*(vertices[i].y - vertices[v_vertices[j]].y) + (vertices[i].z - vertices[v_vertices[j]].z)*(vertices[i].z - vertices[v_vertices[j]].z));
					temp.index = v_vertices[j];
					s_distance.insert(temp);
				}
			}


			vector<distance_vertices> v_distance;
			vector<distance_vertices> temp_distance;
			
			temp_distance.resize(s_distance.size());
			copy(s_distance.begin(), s_distance.end(), temp_distance.begin());
			v_distance.insert(v_distance.begin(), temp_distance.begin(), temp_distance.begin() + (6 - linked_vertices.size()));
			v_vertices_mindistance_vertices.insert(v_vertices_mindistance_vertices.begin(), linked_vertices.begin(), linked_vertices.end());
			for (int k = 0; k < v_distance.size(); k++)
			{
				
				v_vertices_mindistance_vertices.push_back(v_distance[k].index);
			}
		}
		
		if (v_vertices_mindistance_vertices.size() < 4) continue;
		//����������
		//������������ĵ㣬�������������������������ɿռ���
		//�Ȱ�����������ĵ����������γ������ߣ��Ե�һ������ĵ������Ϊ��׼�ߣ��û�׼�߲�������������ٵ�˿ռ�����ȷ���������ڵĲ��棬�õ��ȷ���нǣ��ۺ�����������ȷ�����������ڵ�λ�ã������������
		vector<C3DPoint> temp_triangles(mesh_triangles(i, v_vertices_mindistance_vertices));
		for (int k = 0; k < temp_triangles.size(); k++)
		{
			v_vertices_linked_vertices[temp_triangles[k].x].insert(temp_triangles[k].y);
			v_vertices_linked_vertices[temp_triangles[k].x].insert(temp_triangles[k].z);
			v_vertices_linked_vertices[temp_triangles[k].y].insert(temp_triangles[k].x);
			v_vertices_linked_vertices[temp_triangles[k].y].insert(temp_triangles[k].z);
			v_vertices_linked_vertices[temp_triangles[k].z].insert(temp_triangles[k].x);
			v_vertices_linked_vertices[temp_triangles[k].z].insert(temp_triangles[k].y);

		}
	v_output_triangles.insert(v_output_triangles.end(), temp_triangles.begin(), temp_triangles.end());
	}




	/***********************************************************************************/
	//�����ظ���������
	vector<C3DPoint>v_final_triangles;//���������������������Ƭ
	set<C3DPoint>s_final_triangles;
	for (int n = 0; n < v_output_triangles.size(); n++)
	{
		set<int> s_temp;
		s_temp.insert(v_output_triangles[n].x);
		s_temp.insert(v_output_triangles[n].y);
		s_temp.insert(v_output_triangles[n].z);
		vector<int> v_temp;
		v_temp.resize(s_temp.size());
		copy(s_temp.begin(), s_temp.end(), v_temp.begin());
		C3DPoint p;
		
		p.x = v_temp[0];
		p.y = v_temp[1];
		p.z = v_temp[2];
		s_final_triangles.insert(p);
	}


	v_final_triangles.resize(s_final_triangles.size());
	copy(s_final_triangles.begin(), s_final_triangles.end(), v_final_triangles.begin());




	/*******************************************************/
	//�ڶ����Ż�
	set<line_triangles> s_line;//�����߶�
	vector<line_triangles> v_line;//�����߶�
	vector<line_triangles>   v_line_triangles_combine;//�����Ӧ�߶����ڵ�������
	

	for (int i = 0; i < v_final_triangles.size(); i++)
	{
		line_triangles p;
		line_triangles q;
		line_triangles w;
		p.a = v_final_triangles[i].x;
		p.b = v_final_triangles[i].y;
		p.v_index.push_back(i);
		q.a = v_final_triangles[i].x;
		q.b = v_final_triangles[i].z;
		q.v_index.push_back(i);
		w.a = v_final_triangles[i].y;
		w.b = v_final_triangles[i].z;
		w.v_index.push_back(i);
		s_line.insert(p);
		s_line.insert(q);
		s_line.insert(w);

	}
	v_line.resize(s_line.size());
	copy(s_line.begin(), s_line.end(), v_line.begin());

	line_triangles temp_line;
	temp_line.a = v_line[0].a;
	temp_line.b = v_line[0].b;

	for (int i = 0; i < v_line.size(); i++)
	{
		if (temp_line.a == v_line[i].a&& temp_line.b == v_line[i].b)
		{
			temp_line.v_index.push_back(v_line[i].v_index[0]);
			if (i == (v_line.size() - 1))
				v_line_triangles_combine.push_back(temp_line);
		}
		else if (!(temp_line.a == v_line[i].a&& temp_line.b == v_line[i].b))
		{
			v_line_triangles_combine.push_back(temp_line);
			
			//���ò���
			vector<int>().swap(temp_line.v_index);
			
			temp_line.a = v_line[i].a;
			temp_line.b = v_line[i].b;
			temp_line.v_index.push_back(v_line[i].v_index[0]);
			if (i == (v_line.size() - 1))
				v_line_triangles_combine.push_back(temp_line);
		}
	}

	
	/*******************************************************************************/
	//����Ӧ��������ε��߶���ȡ������ȡ�������߶����������֮����С��ǰ�������㣬ɾ��ʣ���붥������ɵ������Σ����ô˷���
	//����Ӧ��������ε��߶���ȡ������ȡ�������߶�������н�����ǰ�������㣬ɾ��ʣ���붥������ɵ�������



	set<int> s_delete_triangles;
	vector<int> v_delete_triangles;//�����ɾ�������ε�����
	for (int i = 0; i < v_line_triangles_combine.size(); i++)
	{
		if (v_line_triangles_combine[i].v_index.size() >= 3)
		{
			vector<triangle_one_vertice> v_triangle_one_vertice;
			vector<triangle_one_vertice> v_triangle_one_vertice_angle;
			set<triangle_one_vertice> s_triangle_one_vertice;
			vector<triangle_one_vertice> v_sorted_triangle_one_vertice;
			for(int j = 0; j < v_line_triangles_combine[i].v_index.size(); j++)
			{
				triangle_one_vertice temp;
				temp.triangle = v_line_triangles_combine[i].v_index[j];
				temp.vertice = 0;
				temp.angle = 0;
				
				
					if (v_final_triangles[temp.triangle].x != v_line_triangles_combine[i].a&&v_final_triangles[temp.triangle].x != v_line_triangles_combine[i].b)
						temp.vertice = v_final_triangles[temp.triangle].x;
					else if(v_final_triangles[temp.triangle].y != v_line_triangles_combine[i].a&&v_final_triangles[temp.triangle].y != v_line_triangles_combine[i].b)
						temp.vertice = v_final_triangles[temp.triangle].y;
					else
						temp.vertice = v_final_triangles[temp.triangle].z;
				
				v_triangle_one_vertice.push_back(temp);
			}
			v_triangle_one_vertice_angle.resize(v_triangle_one_vertice.size());
			for (int n = 0; n < v_triangle_one_vertice.size(); n++)
			{
				double x1, y1, z1, x2, y2, z2, x3, y3, z3;
				x1 = vertices[v_line_triangles_combine[i].a].x;
				y1 = vertices[v_line_triangles_combine[i].a].y;
				z1 = vertices[v_line_triangles_combine[i].a].z;
				x2 = vertices[v_line_triangles_combine[i].b].x;
				y2 = vertices[v_line_triangles_combine[i].b].y;
				z2 = vertices[v_line_triangles_combine[i].b].z;
				x3 = vertices[v_triangle_one_vertice[n].vertice].x;
				y3 = vertices[v_triangle_one_vertice[n].vertice].y;
				z3 = vertices[v_triangle_one_vertice[n].vertice].z;
				double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
				double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3) + (z1 - z3)*(z1 - z3));
				double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3) + (z2 - z3)*(z2 - z3));
				double temp_angle= (acos((b*b + c*c - a*a) / (2 * b*c)) / 3.1415926) * 180;
				v_triangle_one_vertice_angle[n].angle = temp_angle;
				v_triangle_one_vertice_angle[n].triangle = v_triangle_one_vertice[n].triangle;
				v_triangle_one_vertice_angle[n].vertice = v_triangle_one_vertice[n].vertice;
				s_triangle_one_vertice.insert(v_triangle_one_vertice_angle[n]);
			}
			v_sorted_triangle_one_vertice.resize(s_triangle_one_vertice.size());
			copy(s_triangle_one_vertice.begin(), s_triangle_one_vertice.end(), v_sorted_triangle_one_vertice.begin());
			v_sorted_triangle_one_vertice.erase(v_sorted_triangle_one_vertice.begin(), v_sorted_triangle_one_vertice.begin()+2);
			for (int m = 0; m < v_sorted_triangle_one_vertice.size();m++)
				s_delete_triangles.insert(v_sorted_triangle_one_vertice[m].triangle);
		}
	}
	v_delete_triangles.resize(s_delete_triangles.size());
	copy(s_delete_triangles.begin(), s_delete_triangles.end(), v_delete_triangles.begin());
	for (int i = v_delete_triangles.size() - 1; i >= 0; i--)
	{
		v_final_triangles.erase(v_final_triangles.begin() + v_delete_triangles[i]);
	}
	/**********************************************************************************************/
	//����ֻ��һ�������ε��߶Σ�ȡ�����߶ε��������������ӵĶ��㣬��������һ�����߶�������н�����һ�����㣬�������һ������
	//set<line_triangles> s_line;//�����߶�
	//vector<line_triangles> v_line;//�����߶�
	//vector<line_triangles>   v_line_triangles_combine;//�����Ӧ�߶����ڵ�������

	vector<line_triangles>().swap(v_line);
	vector<line_triangles>().swap(v_line_triangles_combine);
	s_line.clear();
	for (int i = 0; i < v_final_triangles.size(); i++)
	{
		line_triangles p;
		line_triangles q;
		line_triangles w;
		p.a = v_final_triangles[i].x;
		p.b = v_final_triangles[i].y;
		p.v_index.push_back(i);
		q.a = v_final_triangles[i].x;
		q.b = v_final_triangles[i].z;
		q.v_index.push_back(i);
		w.a = v_final_triangles[i].y;
		w.b = v_final_triangles[i].z;
		w.v_index.push_back(i);
		s_line.insert(p);
		s_line.insert(q);
		s_line.insert(w);

	}
	v_line.resize(s_line.size());
	copy(s_line.begin(), s_line.end(), v_line.begin());


	temp_line.a = v_line[0].a;
	temp_line.b = v_line[0].b;
	vector<int>().swap(temp_line.v_index);
	for (int i = 0; i < v_line.size(); i++)
	{
		if (temp_line.a == v_line[i].a&& temp_line.b == v_line[i].b)
		{
			temp_line.v_index.push_back(v_line[i].v_index[0]);
			if (i == (v_line.size() - 1))
				v_line_triangles_combine.push_back(temp_line);
		}
		else if (!(temp_line.a == v_line[i].a&& temp_line.b == v_line[i].b))
		{
			v_line_triangles_combine.push_back(temp_line);

			//���ò���
			vector<int>().swap(temp_line.v_index);

			temp_line.a = v_line[i].a;
			temp_line.b = v_line[i].b;
			temp_line.v_index.push_back(v_line[i].v_index[0]);
			if (i == (v_line.size() - 1))
				v_line_triangles_combine.push_back(temp_line);
		}
	}

	int sum1 = 0;
	int sum2 = 0;
	int sum3 = 0;
	int sum4 = 0;
	for (int i = 0; i < v_line_triangles_combine.size(); i++)
	{
		if (v_line_triangles_combine[i].v_index.size() == 1)
			sum1++;
		if (v_line_triangles_combine[i].v_index.size() == 2)
			sum2++;
		if (v_line_triangles_combine[i].v_index.size() == 3)
			sum3++;
		if (v_line_triangles_combine[i].v_index.size() == 4)
			sum4++;

	}
	//vector<vector<int>>   v_vertices_triangles;//�����Ӧ�±�����ڵ�������
	//vector<set<int>>   v_vertices_linked_vertices;//�������Ӧ�±���Ѿ����ߵĵ������
	vector<vector<int>>().swap(v_vertices_triangles);
	vector<set<int>>().swap(v_vertices_linked_vertices);
	v_vertices_triangles.resize(vertices.size() - 4);

	v_vertices_linked_vertices.resize(vertices.size() - 4);

	for (int i = 0; i < v_final_triangles.size(); i++)
	{
		v_vertices_triangles[v_final_triangles[i].x].push_back(i);
		v_vertices_triangles[v_final_triangles[i].y].push_back(i);
		v_vertices_triangles[v_final_triangles[i].z].push_back(i);

	}

	for (int i = 0; i < v_vertices_triangles.size(); i++)
	{

		set<int> s_vertices;
		vector<int> v_vertices;
		vector<int>linked_vertices;
		if (v_vertices_triangles[i].size() == 0)
		{
			continue;
		}//���ĳ�㲻����κ������Σ���Ϊ������

		for (int j = 0; j < v_vertices_triangles[i].size(); j++)
		{
			s_vertices.insert(v_final_triangles[(v_vertices_triangles[i])[j]].x);
			s_vertices.insert(v_final_triangles[(v_vertices_triangles[i])[j]].y);
			s_vertices.insert(v_final_triangles[(v_vertices_triangles[i])[j]].z);
		}
		s_vertices.erase(i);
		v_vertices.resize(s_vertices.size());
		copy(s_vertices.begin(), s_vertices.end(), v_vertices.begin());
		v_vertices_linked_vertices[i].insert(v_vertices.begin(), v_vertices.end());
	}



	//v_line_triangles_combine.size()

	for (int i = 0; i < v_line_triangles_combine.size() ; i++)
	{
		if (v_line_triangles_combine[i].v_index.size() == 1)
		{
			int triangle = v_line_triangles_combine[i].v_index[0];
			int first_vertice;
			int second_vertice;
			int third_vertice;
			int temp_vertice=0;//�ѹ��������εĵ�
			first_vertice = v_line_triangles_combine[i].a;
			second_vertice = v_line_triangles_combine[i].b;
			
			
				if (v_final_triangles[triangle].x != first_vertice&&v_final_triangles[triangle].x != second_vertice)
					temp_vertice = v_final_triangles[triangle].x;
				else if (v_final_triangles[triangle].y != first_vertice&&v_final_triangles[triangle].y != second_vertice)
					temp_vertice = v_final_triangles[triangle].y;
				else
					temp_vertice = v_final_triangles[triangle].z;
			

			set<int> s_first_vertice_linked_vertices(v_vertices_linked_vertices[first_vertice]);
			s_first_vertice_linked_vertices.erase(second_vertice);
			s_first_vertice_linked_vertices.erase(temp_vertice);
			set<int> s_second_vertice_linked_vertices(v_vertices_linked_vertices[second_vertice]);
			s_second_vertice_linked_vertices.erase(temp_vertice);
			s_second_vertice_linked_vertices.erase(first_vertice);
			//���õ�һ�ַ���ȷ�������㣬����Ѱ�ҵ�һ�������㼯��ڶ��������㼯���غϵ㣬���õ���ڣ����ô˵���Ϊ�ڶ��������εĵ�����

			vector<int> v_first_vertice_linked_vertices;
			v_first_vertice_linked_vertices.resize(s_first_vertice_linked_vertices.size());
			copy(s_first_vertice_linked_vertices.begin(), s_first_vertice_linked_vertices.end(), v_first_vertice_linked_vertices.begin());
			bool temp_bool = FALSE;
			for (int j = 0; j < v_first_vertice_linked_vertices.size(); j++)
			{
				int size = s_second_vertice_linked_vertices.size();
				s_second_vertice_linked_vertices.insert(v_first_vertice_linked_vertices[j]);
				if (size == s_second_vertice_linked_vertices.size())
				{
					third_vertice = v_first_vertice_linked_vertices[j];
					temp_bool = TRUE;
				}
			}
			if (temp_bool == TRUE)
			{
				set<int> ss_triangles;
				vector<int> vv_triangles;
				ss_triangles.insert(first_vertice);
				ss_triangles.insert(second_vertice);
				ss_triangles.insert(third_vertice);
				vv_triangles.resize(ss_triangles.size());
				copy(ss_triangles.begin(), ss_triangles.end(), vv_triangles.begin());
				C3DPoint temp;
				temp.x = vv_triangles[0];
				temp.y = vv_triangles[1];
				temp.z = vv_triangles[2];
				s_final_triangles.insert(temp);
			}
			//�����һ�ַ������У����õڶ��֣����ѵ�һ�������㼯��ڶ��������㼯���һ�����ϣ�Ѱ��ȡ�������߶�������н����Ķ��㣬���������
			else
			{
				vector<int> v_second_vertice_linked_vertices;
				v_second_vertice_linked_vertices.resize(s_second_vertice_linked_vertices.size());
				copy(s_second_vertice_linked_vertices.begin(), s_second_vertice_linked_vertices.end(), v_second_vertice_linked_vertices.begin());
				set<angle_vertices> s_angle_vertices;
				vector<angle_vertices> v_angle_vertices;
				for (int k = 0; k < v_second_vertice_linked_vertices.size(); k++)
				{
					angle_vertices temp;
					temp.index = v_second_vertice_linked_vertices[k];
					double x1, y1, z1, x2, y2, z2, x3, y3, z3;
					x1 = vertices[first_vertice].x;
					y1 = vertices[first_vertice].y;
					z1 = vertices[first_vertice].z;
					x2 = vertices[second_vertice].x;
					y2 = vertices[second_vertice].y;
					z2 = vertices[second_vertice].z;
					x3 = vertices[v_second_vertice_linked_vertices[k]].x;
					y3 = vertices[v_second_vertice_linked_vertices[k]].y;
					z3 = vertices[v_second_vertice_linked_vertices[k]].z;
					double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
					double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3) + (z1 - z3)*(z1 - z3));
					double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3) + (z2 - z3)*(z2 - z3));
					double temp_angle = (acos((b*b + c*c - a*a) / (2 * b*c)) / 3.1415926) * 180;
					temp.angle = temp_angle;

					s_angle_vertices.insert(temp);
				}
				v_angle_vertices.resize(s_angle_vertices.size());
				copy(s_angle_vertices.begin(), s_angle_vertices.end(), v_angle_vertices.begin());
				third_vertice = v_angle_vertices[v_angle_vertices.size() - 1].index;
				set<int> ss_triangles;
				vector<int> vv_triangles;
				ss_triangles.insert(first_vertice);
				ss_triangles.insert(second_vertice);
				ss_triangles.insert(third_vertice);
				vv_triangles.resize(ss_triangles.size());
				copy(ss_triangles.begin(), ss_triangles.end(), vv_triangles.begin());
				C3DPoint temp;
				temp.x = vv_triangles[0];
				temp.y = vv_triangles[1];
				temp.z = vv_triangles[2];
				s_final_triangles.insert(temp);
			}
		}

	}
	vector<C3DPoint>().swap(v_final_triangles);
	v_final_triangles.resize(s_final_triangles.size());
	copy(s_final_triangles.begin(), s_final_triangles.end(), v_final_triangles.begin());
	/*******************************************************************************/
	/**************************************************************************************************/
	set<line_triangles> s_line_2;//�����߶�
	vector<line_triangles> v_line_2;//�����߶�
	vector<line_triangles>   v_line_triangles_combine_2;//�����Ӧ�߶����ڵ�������


	for (int i = 0; i < v_final_triangles.size(); i++)
	{
		line_triangles p;
		line_triangles q;
		line_triangles w;
		p.a = v_final_triangles[i].x;
		p.b = v_final_triangles[i].y;
		p.v_index.push_back(i);
		q.a = v_final_triangles[i].x;
		q.b = v_final_triangles[i].z;
		q.v_index.push_back(i);
		w.a = v_final_triangles[i].y;
		w.b = v_final_triangles[i].z;
		w.v_index.push_back(i);
		s_line_2.insert(p);
		s_line_2.insert(q);
		s_line_2.insert(w);

	}
	v_line_2.resize(s_line_2.size());
	copy(s_line_2.begin(), s_line_2.end(), v_line_2.begin());

	line_triangles temp_line_2;
	temp_line_2.a = v_line_2[0].a;
	temp_line_2.b = v_line_2[0].b;

	for (int i = 0; i < v_line_2.size(); i++)
	{
		if (temp_line_2.a == v_line_2[i].a&& temp_line_2.b == v_line_2[i].b)
		{
			temp_line_2.v_index.push_back(v_line_2[i].v_index[0]);
			if (i == (v_line_2.size() - 1))
				v_line_triangles_combine_2.push_back(temp_line_2);
		}
		else if (!(temp_line_2.a == v_line_2[i].a&& temp_line_2.b == v_line_2[i].b))
		{
			v_line_triangles_combine_2.push_back(temp_line_2);

			//���ò���
			vector<int>().swap(temp_line_2.v_index);

			temp_line_2.a = v_line_2[i].a;
			temp_line_2.b = v_line_2[i].b;
			temp_line_2.v_index.push_back(v_line_2[i].v_index[0]);
			if (i == (v_line_2.size() - 1))
				v_line_triangles_combine_2.push_back(temp_line_2);
		}
	}

/****************************************************************************************************************/
	 sum1 = 0;
	sum2 = 0;
	 sum3 = 0; 
	 sum4 = 0;
	for (int i = 0; i < v_line_triangles_combine_2.size(); i++)
	{
		if (v_line_triangles_combine_2[i].v_index.size() == 1)
			sum1++;
		if (v_line_triangles_combine_2[i].v_index.size() == 2)
			sum2++;
		if (v_line_triangles_combine_2[i].v_index.size() == 3)
			sum3++;
		if (v_line_triangles_combine_2[i].v_index.size() == 4)
			sum4++;

	}




	
	return v_final_triangles;
	//return v_output_triangles;
	 
}