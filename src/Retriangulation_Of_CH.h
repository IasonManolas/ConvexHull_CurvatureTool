//Epanatrigwnopoihsh
Vec3d midOfTwoPoints(Vec3d a, Vec3d b)
{
	Vec3d lineSegAB = subVec(b, a);
	Vec3d mid = a.add(lineSegAB.scale(0.5));
	return mid;
}
void Simple3DScene::CreateFourTriangles(Vec3d a, Vec3d b, Vec3d c, Vec3d mid1, Vec3d mid2, Vec3d mid3, vector<Vec3d> &Tr1, vector<Vec3d>& Tr2, vector<Vec3d> &Tr3, vector<Vec3d>& Tr4)
{
	vector<Vec3d> tempVec;
	//first
	tempVec.push_back(mid1);
	tempVec.push_back(mid2);
	tempVec.push_back(b);
	Tr1 = tempVec;
	//second
	tempVec.clear();
	
	tempVec.push_back(mid2);
	tempVec.push_back(mid3);
	tempVec.push_back(c);
	Tr2 = tempVec;
	//third
	tempVec.clear();

	tempVec.push_back(mid3);
	tempVec.push_back(mid1);
	tempVec.push_back(a);
	Tr3 = tempVec;
	//fourth
	tempVec.clear();

	tempVec.push_back(mid3);
	tempVec.push_back(mid2);
	tempVec.push_back(mid1);
	Tr4 = tempVec;
}

void Simple3DScene::Retriangulation_Of_CH()
{
	static bool FLAG_FIRST_PASS = true;

	if (FLAG_FIRST_PASS)
	{
		vector<Vec3d> &verts = m_model.getVertices();
		for (int i = 0; i < CH.size(); i++)
		{
			Vec3d v1 = verts[CH[i][0]], v2 = verts[CH[i][1]], v3 = verts[CH[i][2]];
			Vec3d mid1 = midOfTwoPoints(v1, v2), mid2 = midOfTwoPoints(v2, v3), mid3 = midOfTwoPoints(v3, v1);
			vector<Vec3d> Tr1, Tr2, Tr3, Tr4;
			CreateFourTriangles(v1, v2, v3, mid1, mid2, mid3, Tr1, Tr2, Tr3, Tr4);
			triangulatedCH.push_back(Tr1);
			triangulatedCH.push_back(Tr2);
			triangulatedCH.push_back(Tr3);
			triangulatedCH.push_back(Tr4);
		}

	}
	FLAG_FIRST_PASS = false;
}
void Simple3DScene::ReRetriangulate()
{
		int N = triangulatedCH.size();
		for (int i = 0; i < N; i++)
		{
			Vec3d v1 = triangulatedCH[i][0], v2 = triangulatedCH[i][1], v3 = triangulatedCH[i][2];
			Vec3d mid1 = midOfTwoPoints(v1, v2), mid2 = midOfTwoPoints(v2, v3), mid3 = midOfTwoPoints(v3, v1);
			vector<Vec3d> Tr1, Tr2, Tr3, Tr4;
			CreateFourTriangles(v1, v2, v3, mid1, mid2, mid3, Tr1, Tr2, Tr3, Tr4);
			triangulatedCH.push_back(Tr1);
			triangulatedCH.push_back(Tr2);
			triangulatedCH.push_back(Tr3);
			triangulatedCH.push_back(Tr4);
		}
		triangulatedCH.erase(triangulatedCH.begin(), triangulatedCH.begin() + N - 1);
}