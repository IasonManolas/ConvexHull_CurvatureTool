Vec3d pointOnLineSeg(Vec3d a, Vec3d b,double c) //dianysma apo to a sto b. AB
{
	Vec3d lineSegAB = subVec(b, a);
	Vec3d point = a.add(lineSegAB.scale(c));
	return point;
}
void Simple3DScene::k_Daktylios(int megethosDaktyliou,double katwfli)
{
	switch (megethosDaktyliou){
	case 1:diwrthwsh_ena_daktyliou(katwfli); break;
	case 2:diwrthwsh_dyo_daktyliou(katwfli); break;
	case 3:diwrthwsh_tria_daktyliou(katwfli); break;
	}
}

void Simple3DScene::diwrthwsh_ena_daktyliou(double katwfli)
{
	vector<Vec3d> &verts = m_model.getVertices();
	int sumCount = 0;
	for (int i = 0; i < extremaPoints.size(); i++) // Gia kathe akrotato
	{
		vertex v = extremaPoints[i];
		double &curvOfV = verKm[v.indexInVerts], goalCurv = verKm[v.indexInVerts] * katwfli;
		Vec3d pointToCorrect = verts[v.indexInVerts];
		int counter = 0;
		do{
			counter++;
			verts.at(v.indexInVerts) = addVec(verts[v.indexInVerts], scaleVec(verNorm[v.indexInVerts], -counter*0.1));
			curvOfPoint(v.indexInVerts);
			
		} while (counter<100 && curvOfV>goalCurv  );
		if (counter == 100) // Me metakinhsh mono twn akrotatwn h kampylothta toy daktyliou den veltiwnetai allo.
		{
			problematicPoints.push_back(verts[v.indexInVerts]);
			sumCount++;
			cout << sumCount << endl;
		}
		
	}
}

void Simple3DScene::diwrthwsh_dyo_daktyliou(double katwfli)
{
	static bool FLAG_FIRST_PASS = true;
	if (FLAG_FIRST_PASS)
	{
	twoRingNeighbours();
	twoRingExtremePoints();
	}
	FLAG_FIRST_PASS = false;
	vector<Vec3d> &verts = m_model.getVertices();
	int sumCount = 0;
	for (int i = 0; i < twoRingExtrPoints.size(); i++) // Gia kathe akrotato
	{
		double &curvOfV = verKm[twoRingExtrPoints[i]], goalCurv = verKm[twoRingExtrPoints[i]] * katwfli;
		Vec3d pointToCorrect = verts[twoRingExtrPoints[i]];
		int counter = 0;
		do{
			counter++;
			verts.at(twoRingExtrPoints[i]) = addVec(verts[twoRingExtrPoints[i]], scaleVec(verNorm[twoRingExtrPoints[i]], -counter*0.1));
			curvOfPoint(twoRingExtrPoints[i]);

		} while (counter<100 && curvOfV>goalCurv);
		if (counter == 100) // Me metakinhsh mono twn akrotatwn h kampylothta toy daktyliou den veltiwnetai allo.
		{
			problematicPoints.push_back(verts[twoRingExtrPoints[i]]);
			sumCount++;
			cout << sumCount << endl;
		}

	}
}
void Simple3DScene::diwrthwsh_tria_daktyliou(double katwfli)
{
	static bool FLAG_FIRST_PASS = true;
	if (FLAG_FIRST_PASS)
	{
		if (twoRingNeig.size()==0) twoRingNeighbours();
		threeRingNeighbours();
		threeRingExtremePoints();
	}
	FLAG_FIRST_PASS = false;
	vector<Vec3d> &verts = m_model.getVertices();
	int sumCount = 0;
	for (int i = 0; i < threeRingExtrPoints.size(); i++) // Gia kathe akrotato
	{
		double &curvOfV = verKm[threeRingExtrPoints[i]], goalCurv = verKm[threeRingExtrPoints[i]] * katwfli;
		Vec3d pointToCorrect = verts[threeRingExtrPoints[i]];
		int counter = 0;
		do{
			counter++;
			verts.at(threeRingExtrPoints[i]) = addVec(verts[threeRingExtrPoints[i]], scaleVec(verNorm[threeRingExtrPoints[i]], -counter*0.1));
			curvOfPoint(threeRingExtrPoints[i]);

		} while (counter<100 && curvOfV>goalCurv);
		if (counter == 100) // Me metakinhsh mono twn akrotatwn h kampylothta toy daktyliou den veltiwnetai allo.
		{
			problematicPoints.push_back(verts[threeRingExtrPoints[i]]);
			sumCount++;
			cout << sumCount << endl;
		}

	}
}
int Simple3DScene::megethos_Daktyliou()
{
	int megethosDaktyliou; 
	do
	{
		cout << "\nDwste megethos dakyliou:";
		cin >> megethosDaktyliou;
	} while (megethosDaktyliou > 3);
	return megethosDaktyliou;

}
double Simple3DScene::katwfli_Kampylothtas()
{
	double katwfli;
	do
	{
		cout << "\nDwste katwfli :";
		cin >> katwfli;
	} while (katwfli < 0 || katwfli>=1);
	return katwfli;
}
void Simple3DScene::curvOfPoint(int index)
{
	for (int j = 0; j < verNeig[index].size(); j++)
		normalsOfvertices(verNeig[index][j]);
	for (int j = 0; j < verNeig[index].size(); j++)
		tiAndKnOfTi(verNeig[index][j]);
	for (int j = 0; j < verNeig[index].size(); j++)
		abcCoefCacl(verNeig[index][j]);
	normalsOfvertices(index);
	tiAndKnOfTi(index);
	abcCoefCacl(index);
	Km(index);
}
void Simple3DScene::Km(int index)
{
	double a = verABC[index][0], c = verABC[index][2];
	double km = (a + c) / 2;
	km = abs(km);
	verKm.at(index)=km;
}
void Simple3DScene::Kg(int index)
{
	vector<Vec3d> &verts = m_model.getVertices();
	vector<vvr::Triangle>& tris = m_model.getTriangles();

	double a = verABC[index][0], b = verABC[index][1], c = verABC[index][2];
		double kg = (a*c - pow(b, 2)) / 4;
		kg = abs(kg);
		verKg.at(index)=kg;


}

void Simple3DScene::abcCoefCacl(int index)
{
	vector<Vec3d> &verts = m_model.getVertices();
	vector<vvr::Triangle>& tris = m_model.getTriangles();
	int i = index;
		vector<double> abc;
		int indexOfMaxKn = distance(knofti[i].begin(), std::max_element(knofti[i].begin(), knofti[i].end())); //max kata metro?
		Vec3d e1 = verTi[i][indexOfMaxKn];
		double a = knofti[i][indexOfMaxKn];
		abc.push_back(a);
		double a11 = 0, a12 = 0, a21 = 0, a22 = 0, a13 = 0, a23 = 0;
		for (int j = 0; j < verNeig[i].size(); j++)
		{
			Vec3d tiTouJGeitona = verTi[i][j];
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
		verABC.at(index)=abc;
	
}

void Simple3DScene::tiAndKnOfTi(int index) //kapoies koryfes yparxoun dyo fores ston verts?
{
	vector<Vec3d> &verts = m_model.getVertices();
	vector<vvr::Triangle>& tris = m_model.getTriangles();
	int i = index;
		vector<Vec3d> ti;
		vector<double> kn;
		for (int j = 0; j < verNeig[i].size(); j++)
		{
			//Calc of ti
			int a = verNeig[i][j];
			Vec3d sub1 = verts[a], sub2 = verts[i];
			Vec3d PiMeionP = subVec(sub1, sub2);
			double DPtouNmePiMeionP = dotProduct(PiMeionP, verNorm[i]);
			Vec3d ginomenoDPmeNormal = scaleVec(verNorm[i], DPtouNmePiMeionP), t = subVec(PiMeionP, ginomenoDPmeNormal);
			Vec3d norT = norVec(t);
			ti.push_back(norT);

			//Calc of kn(ti)
			Vec3d NiMeionN = subVec(verNorm[verNeig[i][j]], verNorm[i]);
			double DPtouPimeionPmeNiMeionN = dotProduct(PiMeionP, NiMeionN), DPtouPimeionPiMetonEautoTou = dotProduct(PiMeionP, PiMeionP);
			double k = -DPtouPimeionPmeNiMeionN / DPtouPimeionPiMetonEautoTou;
			kn.push_back(k);
		}
		verTi.at(index)=ti;
		knofti.at(index)=kn;
}

void Simple3DScene::normalsOfvertices(int index)
{
	vector<Vec3d> &verts = m_model.getVertices();
	vector<vvr::Triangle> &tris = m_model.getTriangles();
	int i = index;
		Vec3d SumWjEpiNormalj = (0, 0, 0);
		vector<Vec3d> centerOfEachTriangle;
		for (int j = 0; j<neigTri[i].size(); j++)
		{
			Vec3d Gj = tris[neigTri[i][j]].getCenter();
			Vec3d GjMeionP = Gj.sub(verts[i]);
			double Wj = 1 / GjMeionP.length();

			Vec3d normaljTrigwnou = tris[neigTri[i][j]].getNormal();
			Vec3d WjEpiNormalj = normaljTrigwnou.scale(Wj);
			SumWjEpiNormalj.add(WjEpiNormalj);
		}
		verNorm.at(index)=(SumWjEpiNormalj.normalize());
	
}
void Simple3DScene::twoRingNeighbours()
{
	vector<Vec3d> &verts = m_model.getVertices();
	twoRingNeig.resize(verts.size());
	for (int i = 0; i < verts.size(); i++)
		for (int j = 0; j < verNeig[i].size(); j++)
		{
			//Oi ena geitones einai kai dyo geitones.
			twoRingNeig.at(i).push_back(verNeig[i][j]);
			for (int k = 0; k < verNeig[verNeig[i][j]].size();k++)
				if (find(twoRingNeig[i].begin(), twoRingNeig[i].end(), verNeig[verNeig[i][j]][k]) == twoRingNeig[i].end()) 
					twoRingNeig.at(i).push_back(verNeig[verNeig[i][j]][k]);
		}
}
void Simple3DScene::threeRingNeighbours()
{
	vector<Vec3d> &verts = m_model.getVertices();
	threeRingNeig.resize(verts.size());
	for (int i = 0; i < verts.size(); i++)
		for (int j = 0; j < twoRingNeig[i].size(); j++)
		{
		//Oi ena geitones einai kai dyo geitones.
		threeRingNeig.at(i).push_back(twoRingNeig[i][j]);
		for (int k = 0; k < verNeig[twoRingNeig[i][j]].size(); k++)
			if (find(threeRingNeig[i].begin(), threeRingNeig[i].end(), verNeig[twoRingNeig[i][j]][k]) == threeRingNeig[i].end())
				threeRingNeig.at(i).push_back(verNeig[twoRingNeig[i][j]][k]);
		}
}
void Simple3DScene::twoRingExtremePoints()
{
	vector<Vec3d> &verts = m_model.getVertices();
	for (int i = 0; i < verts.size(); i++)
		if (isTwoRingLocalExtr(i)) twoRingExtrPoints.push_back(i);
}
bool Simple3DScene::isTwoRingLocalExtr(int index)
{
	double curvOFIndex = verKm[index];
	for (int i = 0; i < twoRingNeig[index].size(); i++) if (verKm[twoRingNeig[index][i]]> curvOFIndex) return false;
	return true;
}
void Simple3DScene::threeRingExtremePoints()
{
	vector<Vec3d> &verts = m_model.getVertices();
	for (int i = 0; i < verts.size(); i++)
		if (isThreeRingLocalExtr(i)) threeRingExtrPoints.push_back(i);
}
bool Simple3DScene::isThreeRingLocalExtr(int index)
{
	double curvOFIndex = verKm[index];
	for (int i = 0; i < threeRingNeig[index].size(); i++) if (verKm[threeRingNeig[index][i]]> curvOFIndex) return false;
	return true;
}
void Simple3DScene::drawTwoRingLocalExtrPoints()
{
	vector<Vec3d> &verts = m_model.getVertices();
	for (int i = 0; i < twoRingExtrPoints.size(); i++)
	{
		Vec3d p(verts[twoRingExtrPoints[i]].x, verts[twoRingExtrPoints[i]].y, verts[twoRingExtrPoints[i]].z);
		Point3D(p.x, p.y, p.z, Colour::blue).draw();
	}
}
void Simple3DScene::drawThreeRingLocalExtrPoints()
{
	vector<Vec3d> &verts = m_model.getVertices();
	for (int i = 0; i < threeRingExtrPoints.size(); i++)
	{
		Vec3d p(verts[threeRingExtrPoints[i]].x, verts[threeRingExtrPoints[i]].y, verts[threeRingExtrPoints[i]].z);
		Point3D(p.x, p.y, p.z, Colour::blue).draw();
	}

}