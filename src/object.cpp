#include "object.h"

// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
	return (int(std::signbit(value)) * (v1-v0)) + v0;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
#include <cmath>

bool Sphere::local_intersect(Ray ray, double t_min, double t_max, Intersection *hit) 
{
    // Sphere centered at (0,0,0) in local space, ray origin transformed to local space.
    double3 oc = center-ray.origin;

    double a = dot(ray.direction, ray.direction);
    double b = -2.0 * dot(ray.direction, oc);
    double c = dot(oc, oc) - radius * radius;

    double discriminant = b*b - 4.0 * a * c ;
    if (discriminant < 0) {
        return false;
    }

    double sqrtd = std::sqrt(discriminant);

    // Find the nearest root within the acceptable range
    double root = (b - sqrtd) / (2.0* a);
	
    if (root <= t_min || t_max <= root) {
        root = (b + sqrtd) / (2.0*a);
        if (root <= t_min || t_max <= root) {
            return false;
        }
    }


    hit->depth = root;
    // Intersection position calculation: position = origin + t * direction
    hit->position = ray.origin + root * ray.direction;
    
    // Calculate the normal vector at the hit position (pointing outward from the sphere surface)
    hit->normal = hit->position / radius;

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	double3 min = world_center - double3{radius, radius, radius};
	double3 max = world_center + double3{radius, radius, radius};
	return AABB{min, max};
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray ray, double t_min, double t_max, Intersection *hit) {
    double denom = dot(normal, ray.direction);
    if (denom > EPSILON){
        double t = dot((center - ray.origin), normal) / denom;

        // abcd
        double3 A = mul(i_transform, double4{})

        hit->normal = normal;
        hit->depth = t;
        hit->position = ray.origin + t * ray.direction;

        return (t >= 0);
    } 
    return false;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
    // Define the quad's vertices in local space (centered around the origin)
    double3 vertices[4] = {
        { half_size,  half_size, 0},
        {-half_size,  half_size, 0},
        {-half_size, -half_size, 0},
        { half_size, -half_size, 0}
    };

    // Transform each vertex to world space
    double3 world_vertices[4];
    for (int i = 0; i < 4; ++i) {
        world_vertices[i] = mul(transform, double4{vertices[i], 1.0}).xyz();
    }

    // Initialize min and max to the first transformed vertex
    double3 min = world_vertices[0];
    double3 max = world_vertices[0];

    // Find the minimum and maximum coordinates
    for (int i = 1; i < 4; ++i) {
        min.x = std::min(min.x, world_vertices[i].x);
        min.y = std::min(min.y, world_vertices[i].y);
        min.z = std::min(min.z, world_vertices[i].z);

        max.x = std::max(max.x, world_vertices[i].x);
        max.y = std::max(max.y, world_vertices[i].y);
        max.z = std::max(max.z, world_vertices[i].z);
    }

    // Return the AABB with the computed min and max bounds
    return AABB{min, max};
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Cylinder::local_intersect(Ray ray, double t_min, double t_max, Intersection *hit) {
    // Cylinder centered at the origin, aligned along the y-axis, with radius and half-height
    double radius_squared = radius * radius;

    // Decompose ray origin and direction for convenience
    double ox = ray.origin.x;
    double oz = ray.origin.z;
    double dx = ray.direction.x;
    double dz = ray.direction.z;
    
    // Intersection with infinite cylinder's side (ignoring y-bounds)
    double a = dx * dx + dz * dz;
    double b = 2 * (ox * dx + oz * dz);
    double c = ox * ox + oz * oz - radius_squared;

    // Solve the quadratic equation a*t^2 + b*t + c = 0
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false; // No real roots; ray misses cylinder

    double sqrt_disc = std::sqrt(discriminant);
    double t1 = (-b - sqrt_disc) / (2 * a);
    double t2 = (-b + sqrt_disc) / (2 * a);

    // Ensure t1 is the smaller root
    if (t1 > t2) std::swap(t1, t2);

    // Check if intersections with the cylinder's side are within y bounds
    double y1 = ray.origin.y + t1 * ray.direction.y;
    if (contains(-half_height, half_height, y1) && contains(t_min, t_max, t1)) {
        hit->depth = t1;
        hit->position = ray.origin + t1 * ray.direction;
        hit->normal = normalize(double3(hit->position.x, 0, hit->position.z));
        return true;
    }

    double y2 = ray.origin.y + t2 * ray.direction.y;
    if (contains(-half_height, half_height, y2) && contains(t_min, t_max, t2)) {
        hit->depth = t2;
        hit->position = ray.origin + t2 * ray.direction;
        hit->normal = normalize(double3(hit->position.x, 0, hit->position.z));
        return true;
    }

    // Intersections with the caps
    double3 cap_normal = double3(0, 1, 0);
    for (int i = -1; i <= 1; i += 2) { // i = -1 for bottom cap, i = 1 for top cap
        double y_cap = i * half_height;
        double t_cap = (y_cap - ray.origin.y) / ray.direction.y;

        if (contains(t_min, t_max, t_cap)) {
            double3 point = ray.origin + t_cap * ray.direction;
            if ((point.x * point.x + point.z * point.z) < radius*radius) {
                hit->depth = t_cap;
                hit->position = point;
                hit->normal = i * cap_normal;
                return true;
            }
        }
    }

    return false; // No valid intersection
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
	return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
//
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
//
bool Mesh::local_intersect(Ray ray,  double t_min, double t_max, Intersection* hit)
{
	return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray  ray, 
							  double t_min, double t_max,
							  Triangle const tri,
							  Intersection *hit)
{
	// Extrait chaque position de sommet des données du maillage.
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A (Pour faciliter les explications)
	double3 const &p1 = positions[tri[1].pi]; // ou Sommet B
	double3 const &p2 = positions[tri[2].pi]; // ou Sommet C

	// Triangle en question. Respectez la convention suivante pour vos variables.
	//
	//     A
	//    / \
	//   /   \
	//  B --> C
	//
	// Respectez la règle de la main droite pour la normale.

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Pour plus de d'informations sur la géométrie, référez-vous à la classe dans object.hpp.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.

	return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	return Object::compute_aabb();
}