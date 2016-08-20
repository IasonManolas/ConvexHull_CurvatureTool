void Simple3DScene::DrawKg()
{
	double max = *max_element(triKg.begin(), triKg.end()), min = *min_element(triKg.begin(), triKg.end()), Diasthma = max - min, Bhma = Diasthma / 5;// Opws kai sto 2D
	vector<vvr::Triangle>& tris = m_model.getTriangles();
	for (int i = 0; i <tris.size(); i++)
	{
		double& Curv = triKg[i];

		if (min <= Curv && Curv < min + Bhma) Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::cyan).draw();
		else if (min + Bhma <= Curv && Curv < min + 2 * Bhma) Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::green).draw();
		else if (min + 2 * Bhma <= Curv && Curv < min + 3 * Bhma) Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::yellow).draw();
		else if (min + 3 * Bhma <= Curv && Curv < min + 4 * Bhma)Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::orange).draw();
		else /*if (Curv <= min + 4 * Bhma && Curv <= max)*/ Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z, //Afhnei kapoia trigwna axrwmatista otan Km
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::red).draw();



	}
}
void Simple3DScene::DrawKm()
{
	double max = *max_element(triKm.begin(), triKm.end()), min = *min_element(triKm.begin(), triKm.end()), Diasthma = max - min, Bhma = Diasthma / 5;// Opws kai sto 2D
	vector<vvr::Triangle>& tris = m_model.getTriangles();
	for (int i = 0; i <tris.size(); i++)
	{
		double& Curv = triKm[i];

		if (min <= Curv && Curv < min + Bhma) Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::cyan).draw();
		else if (min + Bhma <= Curv && Curv < min + 2 * Bhma) Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::green).draw();
		else if (min + 2 * Bhma <= Curv && Curv < min + 3 * Bhma) Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::yellow).draw();
		else if (min + 3 * Bhma <= Curv && Curv < min + 4 * Bhma)Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z,
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::orange).draw();
		else /*if (Curv <= min + 4 * Bhma && Curv <= max)*/ Triangle3D(tris[i].v1().x, tris[i].v1().y, tris[i].v1().z, //Afhnei kapoia trigwna axrwmatista otan Km
			tris[i].v2().x, tris[i].v2().y, tris[i].v2().z,
			tris[i].v3().x, tris[i].v3().y, tris[i].v3().z, Colour::red).draw();



	}
}
void Simple3DScene::CurvOf3DModel() //Ylopoiw thn methodo pou perigrafetai edw: http://www.zju.edu.cn/jzus/2005/A05S1/A05S121.pdf
{
	verNorm.clear(); //exei megethos oses oi koryfes . Kathe grammh i deixnei to normal dianysma sthn vertex[i] .
	verTi.clear(); // Vriskei gia kathe koryfh i ta dianysmata ti ths koryfhs i me olous tous geitones ths ta vazei se vector kai ayton ton vazei sth thesh verTi[i].
	knofti.clear(); // Exei megethos oso to verts. Gia kathe koryfh i vriskei thn kampylothta kn(ti) gia kathe geitona ths i tis vazei se vector kai ayton sth thesh knofti
	verABC.clear(); //Exei gia kathe koryfh i ena vector triwn stoixeiwn tou opoiou to prwto stoixeio einai to kn(tid) to deutero to b kai to trito to c.
	verKg.clear(); //Vector ston opoio to i stoixeio einai h Gaussian kampylothta ths i koryfhs.
	verKm.clear();	//Vector ston opoio to i stoixeio einai h Messh kampylothta ths i koryfhs.
	triKg.clear();
	triKm.clear();
	static bool FLAG_FIRST_PASS = true;
	if (FLAG_FIRST_PASS)
	{
	oneRingNeighbours(); //Vriskei gia kathe koryfh i sto vector verts to index twn geitonwn sto verts kai to vazei sto verNeig[i]
	}
	FLAG_FIRST_PASS = false; 
	normalsOfvertices(); //Vriskei to katheto dianysma se kathe koryfh tou mesh.
	tiAndKnOfTi();
	abcCoefCacl();
	Kg();
	Km();
	curvOfTris();	


}
void Simple3DScene::curvOfTris()
{	
	vector<Vec3d> &verts = m_model.getVertices();
	vector<vvr::Triangle>& tris=m_model.getTriangles();
	
	for(int i=0;i<tris.size();i++)
	{
		triKg.push_back((log(verKg[tris[i].vi1]) + log(verKg[tris[i].vi2]) + log(verKg[tris[i].vi3])) / 3);
		triKm.push_back((log(verKm[tris[i].vi1]) + log(verKm[tris[i].vi2]) + log(verKm[tris[i].vi3])) / 3);
		
	}
}


void Simple3DScene::Kg()
{
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle>& tris = m_model.getTriangles();
		for (int i = 0; i < verts.size(); i++)
		{
			double a = verABC[i][0], b = verABC[i][1], c = verABC[i][2];
			double kg = (a*c - pow(b, 2)) / 4; 
			kg = abs(kg);
			verKg.push_back(kg);
		}

}
void Simple3DScene::Km()
{
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle>& tris = m_model.getTriangles();
		for (int i = 0; i < verts.size(); i++)
		{
			double a = verABC[i][0], c = verABC[i][2];
			double km=(a+c)/2;
			km = abs(km);
			verKm.push_back(km);
		}
}

void Simple3DScene::abcCoefCacl()
{
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle>& tris = m_model.getTriangles();
		for (int i = 0; i < verts.size(); i++)
		{
			vector<double> abc;
			int indexOfMaxKn = distance(knofti[i].begin(), std::max_element(knofti[i].begin(),knofti[i].end())); //max kata metro?
			Vec3d e1 = verTi[i][indexOfMaxKn];
			double a = knofti[i][indexOfMaxKn];
			abc.push_back(a);
			double a11 = 0, a12 = 0, a21 = 0, a22=0,a13 = 0, a23 = 0;
			for (int j = 0; j < verNeig[i].size(); j++)
			{
				Vec3d tiTouJGeitona=verTi[i][j];
				double DPverTimeE1 = dotProduct(tiTouJGeitona, e1);
				if (DPverTimeE1 >= 1) DPverTimeE1 = 1;
				if (DPverTimeE1 <= -1) DPverTimeE1 = -1;
				double thetai = acos(DPverTimeE1);

 				a11 += pow(cos(thetai), 2)*pow(sin(thetai), 2); //a11 mexri a23 swsta
				a12 += cos(thetai)*pow(sin(thetai), 3);
				a21 += cos(thetai)*pow(sin(thetai), 3);
				a22 += pow(sin(thetai), 4);
				a13 += (knofti[i][j] - a*pow(cos(thetai), 2))*cos(thetai)*sin(thetai);
				a23 += (knofti[i][j] - a*pow(cos(thetai), 2))*pow(sin(thetai), 2);
			}
			double b = (a13*a22 - a23*a12) / (a11*a22 - pow(a12, 2)); // b kai c swsta
			double c = (a11*a23 - a12*a13) / (a11*a22 - pow(a12, 2));
			abc.push_back(b);
			abc.push_back(c);
			verABC.push_back(abc);
		}
}

void Simple3DScene::tiAndKnOfTi() //kapoies koryfes yparxoun dyo fores ston verts?
{
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle>& tris = m_model.getTriangles();
		for (int i = 0; i < verts.size(); i++)
		{
			vector<Vec3d> ti;
			vector<double> kn;
			for (int j = 0; j < verNeig[i].size(); j++)
			{
				//Calc of ti
				int a=verNeig[i][j];
				Vec3d sub1 = verts[a], sub2 = verts[i];
				Vec3d PiMeionP = subVec(sub1,sub2);
				double DPtouNmePiMeionP = dotProduct(PiMeionP, verNorm[i]);
				Vec3d ginomenoDPmeNormal = scaleVec(verNorm[i],DPtouNmePiMeionP), t = subVec(PiMeionP,ginomenoDPmeNormal);
				Vec3d norT=norVec(t);
				ti.push_back(norT);

				//Calc of kn(ti)
				Vec3d NiMeionN = subVec(verNorm[verNeig[i][j]], verNorm[i]);
				double DPtouPimeionPmeNiMeionN = dotProduct(PiMeionP, NiMeionN), DPtouPimeionPiMetonEautoTou = dotProduct(PiMeionP, PiMeionP);
				double k = -DPtouPimeionPmeNiMeionN / DPtouPimeionPiMetonEautoTou;
				kn.push_back(k);
			}
			verTi.push_back(ti);
			knofti.push_back(kn);
		}
}

void Simple3DScene::normalsOfvertices()
{
		vector<Vec3d> &verts = m_model.getVertices();
		vector<vvr::Triangle> &tris=m_model.getTriangles();
		for(int i=0; i<verts.size();i++)
		{
			Vec3d SumWjEpiNormalj=(0,0,0);
			vector<Vec3d> centerOfEachTriangle;
			for (int j = 0; j<neigTri[i].size(); j++)
			{
				Vec3d Gj = tris[neigTri[i][j]].getCenter(); 
				Vec3d GjMeionP = subVec(Gj,verts[i]);
				double Wj = 1. / lengthVec(GjMeionP);
				
				Vec3d normaljTrigwnou = tris[neigTri[i][j]].getNormal();
				Vec3d WjEpiNormalj = scaleVec(normaljTrigwnou,Wj);
				SumWjEpiNormalj=addVec(SumWjEpiNormalj, WjEpiNormalj);
			}
			verNorm.push_back(norVec(SumWjEpiNormalj));
		}
}

void Simple3DScene::oneRingNeighbours() // O(n)
{
	vector<vvr::Triangle>& tris=m_model.getTriangles();
	vector<Vec3d> &verts = m_model.getVertices();
	verNeig.resize(verts.size());
	neigTri.resize(verts.size());
	for(int i=0;i<tris.size();i++)
	{
			//Neighbour Vertices
			if (find(verNeig[tris[i].vi1].begin(), verNeig[tris[i].vi1].end(), tris[i].vi2) == verNeig[tris[i].vi1].end()) verNeig[tris[i].vi1].push_back(tris[i].vi2);
			if (find(verNeig[tris[i].vi1].begin(), verNeig[tris[i].vi1].end(), tris[i].vi3) == verNeig[tris[i].vi1].end()) verNeig[tris[i].vi1].push_back(tris[i].vi3);

			if (find(verNeig[tris[i].vi2].begin(), verNeig[tris[i].vi2].end(), tris[i].vi1) == verNeig[tris[i].vi2].end()) verNeig[tris[i].vi2].push_back(tris[i].vi1);
			if (find(verNeig[tris[i].vi2].begin(), verNeig[tris[i].vi2].end(), tris[i].vi3) == verNeig[tris[i].vi2].end()) verNeig[tris[i].vi2].push_back(tris[i].vi3);

			if (find(verNeig[tris[i].vi3].begin(), verNeig[tris[i].vi3].end(), tris[i].vi1) == verNeig[tris[i].vi3].end()) verNeig[tris[i].vi3].push_back(tris[i].vi1);
			if (find(verNeig[tris[i].vi3].begin(), verNeig[tris[i].vi3].end(), tris[i].vi2) == verNeig[tris[i].vi3].end()) verNeig[tris[i].vi3].push_back(tris[i].vi2);

			//Neighbour Triangles
			if (find(neigTri[tris[i].vi1].begin(), neigTri[tris[i].vi1].end(),i) == neigTri[tris[i].vi1].end()) neigTri[tris[i].vi1].push_back(i);


			if (find(neigTri[tris[i].vi2].begin(), neigTri[tris[i].vi2].end(), i) == neigTri[tris[i].vi2].end()) neigTri[tris[i].vi2].push_back(i);
			

			if (find(neigTri[tris[i].vi3].begin(), neigTri[tris[i].vi3].end(), i) == neigTri[tris[i].vi3].end()) neigTri[tris[i].vi3].push_back(i);

	}
}