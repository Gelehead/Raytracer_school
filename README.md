Oscar Lavolet 20242868

# Raytracer
This project has been done in the IFT3355 class, it should support every feature mentionned [here](a24_ift3355_tp2.pdf).

Ideally, the project is expanded to be a raytracing game engine for small games to be developped with a friend, see [this]() or [that]() as ideas that could be implemented in the future


## Notes about the code
_PLEASE READ THIS BEFORE CORRECTING_

It is now common among developpers to use AI to code and the AI often uses commentary on the code to have the developper understand the program. 

However, I tend to also use commentaries to ease the understanding of my code for myself and any other person coming into the project, which could potentially be misinterpreted as blatant AI usage. Please consider this when correcting.


## Documentation for debugging
 - potential problem sources --> `"DOUBT: [...] "`


### doubts so far 
(also good to keep track of progress)


raytracer.cpp
```
// DOUBT: why is it *2* / resolution instead of 1
double delta_x = 2 / scene.resolution[0];
double delta_y = 2 / scene.resolution[1];
// why are u and v correspondances for y and z ??
double3 uVec = {0,1,0};
double3 vVec = {0,0,1};
```

basic.h
```
// random value between (-1, 1)  
static double rand_double_signed() {
	return (double(rand() / double(RAND_MAX) < 0.5 ? -1 : 1)) * double(rand() / double((RAND_MAX)));
}
```

container.cpp
```
// DOUBT: shouldn't it search for aabb boxes too or is it already in the "objects" variable?
if (objects[i]->intersect(ray, t_min, closest_t, &Ihit)){
    if (Ihit.depth < closest_t){
        closest_t = Ihit.depth;
        *hit = Ihit;
        isHit = true;
    }
}
```

Object.cpp
```
//DOUBT: shouldnt there be some calculations to be made ??
hit->normal = world_normal;
```


### advancements and important parts
#### functionnal
 - None

#### Done ?
 - jittery camera
 - Naive intersection
 - sphere collision and AABB definition
 - quad intersection and AABB definition

#### Left to do 
 - Cylinder collision and AABB definition
 - mesh collision and AABB definition
 - AABB 
    - retrieve corners
    - construct aabb
    - compute aabb
    - intersect
 - BVH
 - Shading
    - eclairage local
    - eclairage local 
    - ombres
    - Reflexion
        - reflection mirroir
        - refraction
    - Textures
        - UV mapping
        - Couleur
 - DOF

## Important to understand the code

### Sphere intersection
The formula with "h" coming from the [raytracing in one weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html) :
![simpler sphere intersection](simpler_sphere_intersection.png)




## trucs a faire gaffe 

    shadow acne -> solved by slightly moving normal vector towards the surface
