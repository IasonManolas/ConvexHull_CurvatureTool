void Simple3DScene::Curvature2D()
{
	enterPixelMode();
	vector<double> curv;
	vector<double> lineSegCurv;
	//Ypologismos kampylothtas gia oles tis koryfes toy polygwnou: O(n)
	for (int pi = 0; pi < m_pts.size(); pi++)
	{	

		int a = pi-1;
		if (pi == 0) a = m_pts.size()-1;
		const C2DPoint &p1 = m_pts[a];
		const C2DPoint &p2 = m_pts[pi];
		const C2DPoint &p3 = m_pts[(pi + 1) % m_pts.size()];
		
		C2DTriangle Triangle(p1, p2, p3);
		C2DPoint Centre=Triangle.GetCircumCentre();
		
		LineSeg2D(p1.x, p1.y, p2.x, p2.y,Colour::black).draw();
		
		double R = Centre.Distance(p1);
		curv.push_back(log(1 / R));
	}
	
	
		
		
		//Xrwmatismos kathe koryfhs toy polygwnoy analoga me thn kampylothta tou

		for (int pi = 0; pi < curv.size(); pi++)
		{
			C2DPoint p1 = m_pts[pi], p2 = m_pts[(pi + 1) % m_pts.size()];
			lineSegCurv.push_back((curv[pi] + curv[(pi + 1) % m_pts.size()]) / 2);
		}
		double max = *max_element(lineSegCurv.begin(), lineSegCurv.end()), min = *min_element(lineSegCurv.begin(), lineSegCurv.end()), Diasthma = max - min, Bhma = Diasthma / 5;
		for (int pi = 0; pi < lineSegCurv.size(); pi++)
		{
			C2DPoint p1 = m_pts[pi], p2 = m_pts[(pi + 1) % m_pts.size()];
			if (min <= lineSegCurv[pi] && lineSegCurv[pi] < min + Bhma) LineSeg2D(p1.x, p1.y, p2.x, p2.y, Colour::cyan).draw();
			else if (min + Bhma <= lineSegCurv[pi] && lineSegCurv[pi] < min + 2 * Bhma) LineSeg2D(p1.x, p1.y, p2.x, p2.y, Colour::green).draw();
			else if (min + 2 * Bhma <= lineSegCurv[pi] && lineSegCurv[pi] < min + 3 * Bhma) LineSeg2D(p1.x, p1.y, p2.x, p2.y, Colour::yellow).draw();
			else if (min + 3 * Bhma <= lineSegCurv[pi] && lineSegCurv[pi] < min + 4 * Bhma) LineSeg2D(p1.x, p1.y, p2.x, p2.y, Colour::orange).draw();
			else LineSeg2D(p1.x, p1.y, p2.x, p2.y, Colour::red).draw();
		}
		returnFromPixelMode();
}
void Simple3DScene::drawPolygon2D()
{
    enterPixelMode();

    for (int pi = 0; pi < m_pts.size(); pi++)
    {
        const C2DPoint &p1 = m_pts[pi];
        const C2DPoint &p2 = m_pts[(pi + 1) % m_pts.size()];

        Colour line_col = Colour::black;

        LineSeg2D(p1.x, p1.y, p2.x, p2.y, line_col).draw();
		if (b_show_pts)
			Point2D(p1.x, p1.y, Colour::yellow).draw();
		
    }
	
    returnFromPixelMode();
	
}