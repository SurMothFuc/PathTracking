#pragma once
#include "common.h"
class Scene {
public:
	int  WIDTH =720;
	int HEIGHT =720;
	int  SCREEN_Z = 1;// ͶӰƽ���Z����

	vector<Shape*> shapes;  // ��������ļ���
	

	void addShape(Shape* shape) {
		shapes.push_back(shape);
	}
	
};