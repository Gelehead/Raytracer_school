#pragma once

#include "linalg/linalg.h"
#include <string>
#include <iostream>
using namespace linalg::aliases;


#define PI 3.14159265358979323846
#define EPSILON 1e-6


// ----------------------- vector part ----------------------- 
static double length_squared(double3 u){
	return u.x * u.x + u.y * u.y + u.z * u.z;
}

static double dot(double3 u, double3 v){
	return u.x * v.x + u.y * v.y + u.z + v.z;
}

static double3 le_cross(const double3& a, const double3& b) {
    return double3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// ----------------------- utils -----------------------

static void print(double d){
	std::cout << d << std::endl;
}
static void print(int d){
	std::cout << d << std::endl;
}
static void print(double3 d){
	std::cout << "(" << d.x << ", " << d.y << ", " << d.z << ")" << std::endl;
}
static void print(double4 d){
	std::cout << "(" << d.x << ", " << d.y << ", " << d.z << ", " << d.w << ")" << std::endl;
}
static void print(double4x4 t){
	for (int i = 0; i < 4; i++)
	{
		double4 d = t.row(i);
		print(d);
	}
}
static bool contains(double min, double max, double x){
	return x >= min && x <= max;
}

static bool is_interior(double a, double b) {
    // verifies if both a and b are within the [0,1] range
    if (!contains(0, 1, a) || !contains(0, 1, b))
        return false;
    return true;
}

// ----------------------- random part -----------------------

// Valeur aléatoire entre [0,1)
static double rand_double() {
	return double(rand()) / double((RAND_MAX));
}

// random value between (-1, 1)  
static double rand_double_signed() {
	return (double(rand() / double(RAND_MAX) < 0.5 ? -1 : 1)) * double(rand() / double((RAND_MAX)));
}

// Valeur aléatoire entre [0,1] pour un vecteur
static double2 rand_double2() {
	return double2{rand_double(),rand_double()};
}

static double3 rand_double3(){
	return double3(rand_double(), rand_double(), rand_double());
}

// Valeur aléatoire à l'intérieur d'un disque.
static double2 random_in_unit_disk() {
    while (true) {
        auto p = (2.0 * rand_double2() - 1.0);
        if (length2(p) >= 1) continue;
        return p;
    }
}

// Convertir radian vers degrée
static double rad2deg(double rad) {
	return rad * 360.0 / (2 * PI);
}

// Convertir degrée vers radian
static double deg2rad(double deg) {
	return deg * 2 * PI / 360.0;
}

// Une classe qui représente un rayon
class Ray 
{
public:
	Ray() : origin(0, 0, 0), direction(0, 0, 0) {}
	Ray(double3 origin_, double3 direction_) :
		origin(origin_), direction(direction_)
	{

	}

	double3 origin;    // Origine du rayon
	double3 direction; // Direction du rayon
};