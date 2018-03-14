#ifndef _vertex
#define _vertex
// an array of this struct will hold all vertex information:
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
struct vertex
{
	glm::vec3	position;
	glm::vec3	normal;
	glm::vec3	color;
	glm::vec2	texCoord;
};
#endif