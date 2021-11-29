#pragma once
#include "common.h"
class Scene {
public:
	int  WIDTH =720;
	int HEIGHT =720;
	int  SCREEN_Z = 1;// 投影平面的Z坐标

	vector<Shape*> shapes;  // 几何物体的集合
	

	void addShape(Shape* shape) {
		shapes.push_back(shape);
	}
	
};