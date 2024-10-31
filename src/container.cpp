#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	return true;
}

// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
// 		- Détecter l'intersection avec l'AABB
//			- Si intersection, détecter l'intersection avec la géométrie.
//				- Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	bool isHit = false;
	double closest_t = 0;

	// iterate over every aabb box
	for (int i = 0; i < objects.size(); ++i){
		// if ray collides with aabb
		if (aabbs[i].intersect(ray, t_min, closest_t)){
			// new Intersection object 
			Intersection Ihit;
			// if intersects with any of the objects in the AABB
			// DOUBT: shouldn't it search for aabb boxes too or is it already in the "objects" variable?
			if (objects[i]->intersect(ray, t_min, closest_t, &Ihit)){
				if (Ihit.depth < closest_t){
					closest_t = Ihit.depth;
					*hit = Ihit;
					isHit = true;
				}
			}
		}
	} 
	return isHit;
}