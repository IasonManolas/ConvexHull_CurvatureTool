// Given Data: Data members of instance of Simple3DScene named CH and CHVerts.
void Simple3DScene::curvOfConvexHull()
{
	static bool FLAG_FIRST_PASS = true;
	
	if (FLAG_FIRST_PASS)
	{
		 // vector me megethos CH.size pou deinei gia kathe trigwno tou CH poia shmeia lamvanontai ypopsin gia ton ypologismo ths kampylothtas.
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle>& tris = m_model.getTriangles();
		for (int i = 0; i < CH.size(); i++) // Gia kathe trigwno sto CH vriskw to kontinotero shmeio sto modelo
		{
			math::Triangle tr = CHToMathTr(i);
			vec cenCHtr = tr.Centroid();
			vec cenModelTr = vec3dToVec(tris[0].getCenter());
			int indexMin = 0;
			
			double minDist = cenCHtr.Distance(cenModelTr);
			if (minDist == 0) minDist = 500;
			for (int j = 0; j < tris.size(); j++)
			{
				vec cenModelTr=vec3dToVec(tris[j].getCenter());

				if (cenCHtr.Distance(cenModelTr) < minDist && cenCHtr.Distance(cenModelTr))
				{
					minDist = cenCHtr.Distance(cenModelTr);
					indexMin = j;
				}
			}

			triKmCH.push_back(triKm[indexMin]);
		}
	}

	FLAG_FIRST_PASS = false;
}
void Simple3DScene::DrawKmCH()
{
	double max = *max_element(triKmCH.begin(), triKmCH.end()), min = *min_element(triKmCH.begin(), triKmCH.end()), Diasthma = max - min, Bhma = Diasthma /5; // Opws kai sto 2D
	vector<Vec3d> &verts = m_model.getVertices();
	for (int i = 0; i <CH.size(); i++)
	{
		double& Curv = triKmCH[i];

		Vec3d v1, v2, v3;
		v1 =  verts[CH[i][0]];
		v2 =  verts[CH[i][1]];
		v3 =  verts[CH[i][2]];

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
