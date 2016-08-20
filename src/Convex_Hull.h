void rotateTri(int Tri[3])
{
	int Temp = Tri[1];
	Tri[1] = Tri[2];
	Tri[2] = Temp;
}

bool areAdjacent(vector<int> a, vector<int> b, int &indexAdj1, int &indexAdj2)
{

	if (a[0] == b[0] || a[0] == b[1] || a[0] == b[2])
	{
		if (a[1] == b[0] || a[1] == b[1] || a[1] == b[2])
		{
			indexAdj1 = a[0];
			indexAdj2 = a[1];
			return true;
		}
		else if (a[2] == b[0] || a[2] == b[1] || a[2] == b[2])
		{
			indexAdj1 = a[0];
			indexAdj2 = a[2];
			return true;
		}
	}
	else if (a[1] == b[0] || a[1] == b[1] || a[1] == b[2])
	{
		if (a[0] == b[0] || a[0] == b[1] || a[0] == b[2])
		{
			indexAdj1 = a[0];
			indexAdj2 = a[1];
			return true;
		}
		else if (a[2] == b[0] || a[2] == b[1] || a[2] == b[2])

		{
			indexAdj1 = a[1];
			indexAdj2 = a[2];
			return true;
		}
	}

	return false;

}
void Simple3DScene::getHorizonEdges(vector<vector<int>> adjEdges, int triIndex, vector<vector<int>> &horizonEdges)
{
	if (adjEdges.size() == 2) // one tri edge belongs to horizon
	{
		vector<int> e1;
		int freeVertex1,freeVertex2;

		if ((adjEdges[0][0] == CH[triIndex][0] || adjEdges[0][1] == CH[triIndex][0]) && (adjEdges[1][0] == CH[triIndex][0] || adjEdges[1][1] == CH[triIndex][0]))
		{

			freeVertex1 = CH[triIndex][1];
			freeVertex2 = CH[triIndex][2];
		}
		if ((adjEdges[0][0] == CH[triIndex][1] || adjEdges[0][1] == CH[triIndex][1]) && (adjEdges[1][0] == CH[triIndex][1] || adjEdges[1][1] == CH[triIndex][1]))
		{
			freeVertex1 = CH[triIndex][0];
			freeVertex2 = CH[triIndex][2];
		}
		if ((adjEdges[0][0] == CH[triIndex][2] || adjEdges[0][1] == CH[triIndex][2]) && (adjEdges[1][0] == CH[triIndex][2] || adjEdges[1][1] == CH[triIndex][2]))
		{
			freeVertex1 = CH[triIndex][0];
			freeVertex2 = CH[triIndex][1];
		}
		e1.push_back(freeVertex1);
		e1.push_back(freeVertex2);
		horizonEdges.push_back(e1);
	}
	else if (adjEdges.size() == 1) // two tri edges belong to the horizon
	{
		vector<int> e1, e2;
		int freeVertexIndex= CH[triIndex][0]; //estw
		if (adjEdges[0][0] == freeVertexIndex || adjEdges[0][1] == freeVertexIndex)
			freeVertexIndex = CH[triIndex][1];
		if (adjEdges[0][0] == freeVertexIndex || adjEdges[0][1] == freeVertexIndex)
			freeVertexIndex = CH[triIndex][2];

		e1.push_back(freeVertexIndex);
		e2.push_back(freeVertexIndex);
		e1.push_back(adjEdges[0][0]);
		e2.push_back(adjEdges[0][1]);
		horizonEdges.push_back(e1);
		horizonEdges.push_back(e2);
	}
	else // all tri edges belong to horizon
	{
		vector<int> e1, e2, e3;
		e1.push_back(CH[triIndex][0]);
		e1.push_back(CH[triIndex][1]);
		e2.push_back(CH[triIndex][1]);
		e2.push_back(CH[triIndex][2]);
		e3.push_back(CH[triIndex][2]);
		e3.push_back(CH[triIndex][0]);
		horizonEdges.push_back(e1); 
		horizonEdges.push_back(e2);
		horizonEdges.push_back(e3);
	}
}
void Simple3DScene::ConvexHull3D()
{
	vector<Vec3d> &verts = m_model.getVertices();
	int N = verts.size();
	static bool FLAG_FIRST_PASSa = true;
	static bool FLAG_FIRST_PASS = true;
	if (FLAG_FIRST_PASSa)
	{
	if (FLAG_FIRST_PASS)
	{
	CHInitialization();
	}
	FLAG_FIRST_PASS = false;
	for (int i = 0; i < P.size(); i++)
	{
		CHConflicts(i);
		CHExpand(i);
	}
	}
	FLAG_FIRST_PASSa = false;
	
}

void Simple3DScene::CHConflicts(int i)
{
	vector<Vec3d> &verts = m_model.getVertices();

		Fconf.resize(P.size());

			vec vert(verts[P[i]].x, verts[P[i]].y, verts[P[i]].z);
			for (int j = 0; j < CH.size(); j++) // for all faces of the CH
			{
				math::Triangle Tr = CHToMathTr(j);
				Plane pl=Tr.PlaneCCW();
				if (pl.IsOnPositiveSide(vert) || Tr.PlaneCCW().Distance(vert)<0.000001)  //vert can see j triangle of CH
				{
					Fconf[i].push_back(j); // [p1[0,1,3],p2[2,3],..]
					
				}
			}
}
void Simple3DScene::CHExpand(int indexInPofNewPoint)
{
	vector<Vec3d> &verts = m_model.getVertices();

	
	if (!Fconf[indexInPofNewPoint].empty())
	{
		vector<vector<int>> horizonEdges;
		PiHorizon(indexInPofNewPoint, horizonEdges);
		DeleteVisibleFaces(indexInPofNewPoint);
		CreateNewFaces(indexInPofNewPoint, horizonEdges);
		
	}
	
}
void Simple3DScene::PiHorizon(int i, vector<vector<int>> &horizonEdges) //vector me ta index toy orizonta
{

	for (int j = 0; j < Fconf[i].size(); j++)
	{
		vector<vector<int>> adjEdges;
		
		for (int k = 0; k < Fconf[i].size(); k++)
		{
			if (k == j) continue;
			int indexAdj1, indexAdj2;
			if (areAdjacent(CH[Fconf[i][j]], CH[Fconf[i][k]], indexAdj1, indexAdj2))
			{
				vector<int> temp;
				temp.push_back(indexAdj1);
				temp.push_back(indexAdj2);
				adjEdges.push_back(temp);

				if (adjEdges.size() == 3) break;
			}
		}
		if (adjEdges.size() <3) getHorizonEdges(adjEdges, Fconf[i][j], horizonEdges);
	}

}
void Simple3DScene::DrawHorizon(vector<vector<int>> horizonEdges)
{
	vector<Vec3d> &verts = m_model.getVertices();
	int NumHorizonEdges = horizonEdges.size();
	for (int i = 0; i < NumHorizonEdges; i++)
	{
		LineSeg3D(verts[horizonEdges[i][0]].x, verts[horizonEdges[i][0]].y, verts[horizonEdges[i][0]].z,
			verts[horizonEdges[i][1]].x, verts[horizonEdges[i][1]].y, verts[horizonEdges[i][1]].z, Colour::red).draw();
		
	}
}
void Simple3DScene::DeleteVisibleFaces(int indexOfNewPoint)
{
	sort(Fconf[indexOfNewPoint].begin(), Fconf[indexOfNewPoint].end());
	reverse(Fconf[indexOfNewPoint].begin(), Fconf[indexOfNewPoint].end());

	for (int i = 0; i < Fconf[indexOfNewPoint].size(); i++)
		CH.erase(CH.begin() + Fconf[indexOfNewPoint][i]);

}
void Simple3DScene::CreateNewFaces(int indexOfNewPoint, vector<vector<int>> horizonEdges)
{

	CHVerts.push_back(P[indexOfNewPoint]);
	newCM(P[indexOfNewPoint]);
	for (int i = 0; i < horizonEdges.size(); i++)
	{
		vector<int> tri;
		int Tri[3] = { P[indexOfNewPoint], horizonEdges[i][0], horizonEdges[i][1] };
		
		if (!IsCCW(Tri)) rotateTri(Tri);
		
		tri.assign(Tri, Tri + 3);
		
		CH.push_back(tri);
		
	}

}
void Simple3DScene::DrawCH()
{
	vector<Vec3d> &verts = m_model.getVertices();
	for (int i = 0; i < CH.size(); i++)
		Triangle3D(verts[CH[i][0]].x, verts[CH[i][0]].y, verts[CH[i][0]].z,
		verts[CH[i][1]].x, verts[CH[i][1]].y, verts[CH[i][1]].z,
		verts[CH[i][2]].x, verts[CH[i][2]].y, verts[CH[i][2]].z,Colour::yellow).draw();
	

	
}
void Simple3DScene::DrawConflictsOf(int indexInP)
{
	vector<Vec3d> &verts = m_model.getVertices();

	Point3D(verts[P[indexInP]].x, verts[P[indexInP]].y, verts[P[indexInP]].z, Colour::orange).draw();
	
	for (int i = 0; i < Fconf[indexInP].size(); i++)
		Triangle3D(verts[CH[Fconf[indexInP][i]][0]].x, verts[CH[Fconf[indexInP][i]][0]].y, verts[CH[Fconf[indexInP][i]][0]].z,
		verts[CH[Fconf[indexInP][i]][1]].x, verts[CH[Fconf[indexInP][i]][1]].y, verts[CH[Fconf[indexInP][i]][1]].z,
		verts[CH[Fconf[indexInP][i]][2]].x, verts[CH[Fconf[indexInP][i]][2]].y, verts[CH[Fconf[indexInP][i]][2]].z, Colour::white).draw();
}
void Simple3DScene::CHInitialization() 
{
	static bool FLAG_FIRST_PASS = true;
	vector<Vec3d> &verts = m_model.getVertices();

	static int i1, i2, i3, i4;
	if (FLAG_FIRST_PASS)
	{
		vec p1, p2, p3, p4;
		for (int i = 0; i < verts.size(); i++) P.push_back(i);

		p1.Set(verts[P[0]].x, verts[P[0]].y, verts[P[0]].z);
		i1 = P[0];
		p2.Set(verts[P[1]].x, verts[P[1]].y, verts[P[1]].z);
		i2 = P[1];
		P.erase(P.begin());
		P.erase(P.begin());
		LineSegment lineSeg(p1, p2);
		Line line(lineSeg);
		std::random_shuffle(P.begin(), P.end());
		for (int i = 0; i < P.size(); i++)
		{
			if (i1 == P[i] || i2 == P[i]) continue;
			p3.Set(verts[P[i]].x, verts[P[i]].y, verts[P[i]].z);
			
			if (!line.Contains(p3))
			{
				i3 = P[i];
				Plane plane(p1, p2, p3);
				P.erase(P.begin() + i);
				std::random_shuffle(P.begin(), P.end());
				for (int j = 0; j < P.size(); j++)
				{
					if (P[j] == i1 || P[j] == i2 || P[j] == i3) continue;
					p4.Set(verts[P[j]].x, verts[P[j]].y, verts[P[j]].z);
					
					if (!plane.Distance(p4) == 0)
					{
						i4 = P[j];
						P.erase(P.begin() + j);
						break;
					}
				}
			}
			break;
		}
	// Random Permutation of P
	std::random_shuffle(P.begin(), P.end());
	//found 4 appropriate points
	

	CHVerts.push_back(i1);
	CHCM= verts[i1];

	CHVerts.push_back(i2);
	newCM(i2);

	CHVerts.push_back(i3);
	newCM(i3);

	CHVerts.push_back(i4);
	newCM(i4);


	vector<int> tri1, tri2, tri3, tri4;
	int Tri1[3] = { i1, i2, i3 };
	if (!IsCCW(Tri1)) rotateTri(Tri1);
	tri1.assign(Tri1, Tri1 + 3);
	int Tri2[] = { i4, i3, i2 };
	if (!IsCCW(Tri2)) rotateTri(Tri2);
	tri2.assign(Tri2, Tri2 + 3);
	int Tri3[] = { i4, i1, i3 };
	if (!IsCCW(Tri3)) rotateTri(Tri3);
	tri3.assign(Tri3, Tri3 + 3);
	int Tri4[] = { i1, i2, i4 };
	if (!IsCCW(Tri4)) rotateTri(Tri4);
	tri4.assign(Tri4, Tri4 + 3);
	CH.push_back(tri1);
	CH.push_back(tri2);
	CH.push_back(tri3);
	CH.push_back(tri4);
}
	FLAG_FIRST_PASS = false;
	

}
void Simple3DScene::CHInitialConflicts() // Ο(n) 
{	
	vector<Vec3d> &verts = m_model.getVertices();
	static bool FLAG_FIRST_PASS = true;

	if (FLAG_FIRST_PASS)
	{
	Fconf.resize(P.size());
	Pconf.resize(CH.size());
	
	for (int i = 0; i < P.size(); i++) //for all remaining points 
	{
		vec vert(verts[P[i]].x, verts[P[i]].y, verts[P[i]].z);
		for (int j = 0; j < CH.size(); j++) // for all faces of the CH
		{
			math::Triangle Tr = CHToMathTr(j);
			if (Tr.PlaneCCW().IsOnPositiveSide(vert)) //vert can see j triangle of CH
			{
				Fconf[i].push_back(j); // [p1[0,1,3],p2[2,3],..]
			}
		}
	}
	}
	FLAG_FIRST_PASS = false;
}
void Simple3DScene::newCM(int indexOfNewPoint) //index in verts[]
{
	vector<Vec3d> &verts = m_model.getVertices();
	int N = CHVerts.size();
	CHCM.scale(N-1);
	Vec3d NewPoint = verts[indexOfNewPoint];
	CHCM.add(NewPoint);
	CHCM.scale(1./N);

}
