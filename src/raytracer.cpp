#include "raytracer.h"
#include "basic.h"
#include "container.h"

void Raytracer::render(const Scene& scene, Frame* output)
{       
    // Crée le z_buffer.
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far; //Anciennement DBL_MAX. À remplacer avec la valeur de scene.camera.z_far
    }


	//---------------------------------------------------------------------------------------------------------------
	// Nous vous fournissons ci-bas du code pour une caméra orthographique très simple. Cette caméra peut être utilisée pour tester l’intersection avec la sphère.
	// Vous devez utiliser la scène de test portho.ray pour utiliser cette caméra. 
	// Notez que votre code de caméra ne doit pas être basé sur ce code de caméra. Ce code n’est là que pour prendre en compte le développement initial du test d’intersection.
	// Pour utiliser cette caméra, vous devez supprimer les commentaires qui rendent inactive cette partie du code, et mettre en commentaires la boucle d’image originale.

/* 	CameraOrthographic camOrth;
	double3 uVec{ 0,1,0 };
	double3 vVec{ 0,0,1 };
	double y_shift = 2.0 / scene.resolution[1];
	double x_shift = 2.0 / scene.resolution[0];

	for (int y = 0; y < scene.resolution[1]; y++) {
		if (y % 40) {
			std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
		}

		for (int x = 0; x < scene.resolution[0]; x++) {
			double3 color{ 0,0,0 };

			Intersection hit;

			double3 offset = {-0.3,-0.5,0.5};

			double3 rayOrigin = camOrth.minPosition
			+ uVec * x_shift * x 
			+ vVec * y_shift * y
			+ offset;
			double3 rayDirection{ 1,0,0 };
			Ray ray = Ray(rayOrigin, rayDirection);
			double itHits = 0;

			double z_depth = scene.camera.z_far;
			if (scene.container->intersect(ray, EPSILON, z_depth, &hit)) {
				Material& material = ResourceManager::Instance()->materials[hit.key_material];
				color = hit.normal;
				itHits = 1.0f;
			}

			output->set_color_pixel(x, y, color);
			output->set_depth_pixel(x, y, itHits);
		}
	} */

// ----------------------------------------------------- //

	Camera cam;

	double3 lookAt = normalize(scene.camera.center - scene.camera.position);
	double3 rCcam = normalize(cross(lookAt, scene.camera.up)); 
    double3 uCam = normalize(cross(rCcam, lookAt));

	double3 image_center = lookAt * scene.camera.z_near;

	double image_height = tan(deg2rad(scene.camera.fovy)) * scene.camera.z_near;
	double image_length = image_height / scene.camera.aspect;

	double delta_u = image_length / scene.resolution[0];
	double delta_v = image_height / scene.resolution[1];


	double3 minPos = 
		scene.camera.position
	  + scene.camera.z_near * lookAt
	  - uCam * ((0.5 * image_height) + delta_v) 
	  - rCcam * ((0.5 * image_length) + delta_u)
	;

	// loop over every pixel
	for (int y = 0; y < scene.resolution[1]; y++) {
		if (y % 40) {
			std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
		}

		for (int x = 0; x < scene.resolution[0]; x++) {
			double3 average_color{ 0,0,0 };
			double average_z_depth = 0;

			Intersection hit;

			for (size_t i = 0; i < scene.samples_per_pixel; i++)
			{
				double3 color = {0,0,0};
				double depth = 0;

				double jitterU = (rand_double_signed() / 2) * delta_u;
				double jitterV = (rand_double_signed() / 2) * delta_v;

				double3 jitter = {0,jitterU, jitterV};

				double3 rayOrigin = cam.position;

				double3 rayDirection = normalize(
					minPos 
					+ x * delta_u * rCcam
					+ y * delta_v * uCam
					+ jitter
					- scene.camera.position
				);

				Ray ray = Ray(rayOrigin, rayDirection);

				trace(scene, ray, 0, &color, &depth);

				average_color += color;
				average_z_depth += depth;
			}

			average_color = average_color / scene.samples_per_pixel;
			average_z_depth = average_z_depth / scene.samples_per_pixel;

			if(average_z_depth >= scene.camera.z_near && average_z_depth <= scene.camera.z_far && 
				average_z_depth < z_buffer[x + y*scene.resolution[0]]) {
				z_buffer[x + y*scene.resolution[0]] = average_z_depth;

				std::cout << "hit resolved" << std::endl;

				// Met à jour la couleur de l'image (et sa profondeur)
				output->set_color_pixel(x, y, average_color);
				output->set_depth_pixel(x, y, (average_z_depth - scene.camera.z_near) / 
										(scene.camera.z_far-scene.camera.z_near));
			}
		}
	}
    delete[] z_buffer;

}
// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		- Détermine si le rayon intersecte la géométrie.
//      	- Calculer la contribution associée à la réflexion.
//			- Calculer la contribution associée à la réfraction.
//			- Mettre à jour la couleur avec le shading + 
//			  Ajouter réflexion selon material.reflection +
//			  Ajouter réfraction selon material.refraction 
//            pour la couleur de sortie.
//          - Mettre à jour la nouvelle profondeure.
void Raytracer::trace(const Scene& scene,
					  Ray ray, int ray_depth,
					  double3* out_color, double* out_z_depth)
{
	Intersection hit;
	// Fait appel à l'un des containers spécifiées.
	if(scene.container->intersect(ray, EPSILON, *out_z_depth, &hit)) {

		Material& material = ResourceManager::Instance()->materials[hit.key_material];
		double3 ambiant_color = {0.1, 0.1, 0.1};

		
		//*out_color += material.color_albedo;
		*out_color += hit.normal;
		*out_z_depth += hit.depth; 


		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réflection d'un rayon de manière récursive.
		
		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
		// 
		// Assumez que l'extérieur/l'air a un indice de réfraction de 1.
		//
		// Toutes les géométries sont des surfaces et non pas de volumes.

		// *out_color = 
		// *out_z_depth =
	}
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		* Calculer la contribution des lumières dans la scène.
//			- Itérer sur toutes les lumières.
//				- Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//	          	- Inclure la contribution diffuse. (Faites attention au produit scalare. >= 0)
//   	  	- Inclure la contribution ambiante
//      * Calculer si le point est dans l'ombre
//			- Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occludé.
//				- Ne pas considérer les points plus loins que la lumière.
//			- Par la suite, intégrer la pénombre dans votre calcul
//		* Déterminer la couleur du point d'intersection.
//        	- Si texture est présente, prende la couleur à la coordonnées uv
//			- Si aucune texture, prendre la couleur associé au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit)
{
	// Material& material = ResourceManager::Instance()->materials[hit.key_material]; lorsque vous serez rendu à la partie texture.
	return double3{0,0,0};
}
