#pragma once
class Node
{
public:
	
	long double x;
	long double y;
	long double t;
	bool onBound = false;

	Node(long double x, long double y, bool onBound);
	~Node();

};
