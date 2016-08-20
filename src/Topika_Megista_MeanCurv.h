bool compareByCurv(const vertex &a, const vertex &b)
{
	return a.verKmCurv < b.verKmCurv;
}
void Simple3DScene::Topika_Megista()
{
	vector<Vec3d> &verts = m_model.getVertices();
	static bool FLAG_FIRST_PASS = true;

	if (FLAG_FIRST_PASS)
	{
		for (int i = 0; i < verts.size(); i++)
		{
			vertex v;
			v.indexInVerts = i;
			v.verKmCurv = verKm[i];
			vertsIndexAndCurv.push_back(v);
		}
		sort(vertsIndexAndCurv.begin(), vertsIndexAndCurv.end(), compareByCurv);
		for (int i = 0; i < vertsIndexAndCurv.size(); i++) if (IsLocalExtrema(vertsIndexAndCurv[i])) extremaPoints.push_back(vertsIndexAndCurv[i]);

	}
	FLAG_FIRST_PASS = false;

}
bool Simple3DScene::IsLocalExtrema(vertex vertToCheck)
{
	double KmCurvOfCheckPoint = vertToCheck.verKmCurv;
	for (int i = 0; i < verNeig[vertToCheck.indexInVerts].size(); i++) if (KmCurvOfCheckPoint < verKm[verNeig[vertToCheck.indexInVerts][i]]) return false; // Otan vrei oti kapoios apo tous ena geitones tou shmeiou vertToCheck exei megalyterh kampylothta apo to idio to shmeio epistrefei false. Dhladh to shmeio elegxoy den einai topiko megisto.
	return true;
}
void Simple3DScene::DrawExtremaPoints(int kDaktylios)
{
	vector<Vec3d> &verts = m_model.getVertices();
	if (kDaktylios == 1)
		for (int i = 0; i < extremaPoints.size(); i++)
		{
		Vec3d v = verts[extremaPoints[i].indexInVerts];
		Point3D(v.x, v.y, v.z, Colour::blue).draw();
		}
	else if (kDaktylios == 2)
	{
		static bool FLAG_FIRST_PASS = true;
		if (FLAG_FIRST_PASS)
		{
			twoRingNeighbours();
			twoRingExtremePoints();
		}
		FLAG_FIRST_PASS = false;

	
		for (int i = 0; i < twoRingExtrPoints.size(); i++)
		{
		Vec3d v = verts[twoRingExtrPoints[i]];
		Point3D(v.x, v.y, v.z, Colour::blue).draw();
		}
	}
	else
	{
		static bool FLAG_FIRST_PASS = true;
		if (FLAG_FIRST_PASS)
		{
			if (twoRingNeig.size() == 0) twoRingNeighbours();
			threeRingNeighbours();
			threeRingExtremePoints();
		}
		FLAG_FIRST_PASS = false;
		for (int i = 0; i < threeRingExtrPoints.size(); i++)
		{
			Vec3d v = verts[threeRingExtrPoints[i]];
			Point3D(v.x, v.y, v.z, Colour::blue).draw();
		}
	}
}