//Curvature of RetriangulatedCH
math::Triangle Simple3DScene::triCHtoMathTr(int indexInTriCH)
{
	vector<Vec3d> vec = triangulatedCH[indexInTriCH];
	math::Triangle Tri;
	Tri.a.Set(vec[0].x, vec[0].y, vec[0].z);
	Tri.b.Set(vec[1].x, vec[1].y, vec[1].z);
	Tri.c.Set(vec[2].x, vec[2].y, vec[2].z);
	return Tri;
}
Vec3d centroid(Vec3d a, Vec3d b, Vec3d c)
{
	Vec3d temp = addVec(a, b), temp1 = addVec(temp, c);
	temp1=scaleVec(temp1, 1. / 3.);
	return temp1;
}
void Simple3DScene::Curv_Of_Retri_CH()
{


		reTriKmCH.clear();
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle>& tris = m_model.getTriangles();
		for (int i = 0; i < triangulatedCH.size(); i++) // Gia kathe trigwno sto CH vriskw to kontinotero shmeio sto modelo
		{
			Vec3d v1 = triangulatedCH[i][0], v2 = triangulatedCH[i][1], v3 = triangulatedCH[i][2];

			vec cenCHtr = vec3dToVec(centroid(v1, v2, v3));
			vec cenModelTr = vec3dToVec(tris[0].getCenter());
			int indexMin = 0;

			double minDist = cenCHtr.Distance(cenModelTr);
			if (minDist == 0) minDist = 500;
			for (int j = 0; j < tris.size(); j++)
			{
				vec cenModelTr = vec3dToVec(tris[j].getCenter());

				if (cenCHtr.Distance(cenModelTr) < minDist && cenCHtr.Distance(cenModelTr))
				{
					minDist = cenCHtr.Distance(cenModelTr);
					indexMin = j;
				}
			}

			reTriKmCH.push_back(triKm[indexMin]);
		}

}

void Simple3DScene::DrawRetriKmCH()
{
	double max = *max_element(reTriKmCH.begin(), reTriKmCH.end()), min = *min_element(reTriKmCH.begin(), reTriKmCH.end()), Diasthma = max - min, Bhma = Diasthma / 5; // Opws kai sto 2D
	for (int i = 0; i <triangulatedCH.size(); i++)
	{
		double& Curv = reTriKmCH[i];

		Vec3d v1, v2, v3;
		v1 = triangulatedCH[i][0];
		v2 = triangulatedCH[i][1];
		v3 = triangulatedCH[i][2];
		
		if (min <= Curv && Curv < min + Bhma) Triangle3D(v1.x, v1.y, v1.z,
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z, Colour::cyan).draw();
		else if (min + Bhma <= Curv && Curv < min + 2 * Bhma) Triangle3D(v1.x, v1.y, v1.z,
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z, Colour::green).draw();
		else if (min + 2 * Bhma <= Curv && Curv < min + 3 * Bhma) Triangle3D(v1.x, v1.y, v1.z,
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z, Colour::yellow).draw();
		else if (min + 3 * Bhma <= Curv && Curv < min + 4 * Bhma) Triangle3D(v1.x, v1.y, v1.z,
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z, Colour::orange).draw();
		else /*if (Curv <= min + 4 * Bhma && Curv <= max)*/ Triangle3D(v1.x, v1.y, v1.z,
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z, Colour::red).draw();



	}
}
