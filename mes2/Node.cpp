#include "Node.h"
#include "Input.h"

Node::Node(long double x, long double y, bool onBound)
{
	this->x = x;
	this->y = y;
	this->t = Input::getInitialTemperature();
	this->onBound = onBound;
}


Node::~Node()
{
}

