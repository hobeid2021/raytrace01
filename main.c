#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#define IMG_NAME "out.ppm"
#define IMG_W 500
#define IMG_H 500
#define T_MIN 1
#define T_MAX 100
#define RECURSION_LIMIT 3
#define SPHERE_NUM 4
#define LIGHT_NUM 3

typedef struct {
	float x;
	float y;
	float z;
} vec3;

typedef struct {
	float x;
	float y;
} vec2;

vec3 Add(vec3 a, vec3 b) {
	vec3 v = {a.x + b.x, a.y + b.y, a.z + b.z};
	return v;
}

vec3 Sub(vec3 a, vec3 b) {
	vec3 v = {a.x - b.x, a.y - b.y, a.z - b.z};
	return v;
}

float Dot(vec3 a, vec3 b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

float Mag(vec3 v) {
	return sqrt(Dot(v,v));
}

vec3 AddScalar(vec3 v, float x) {
	return (vec3){v.x+x, v.y+x, v.z+x};
}

vec3 MultScalar(vec3 v, float x) {
	return (vec3){v.x*x, v.y*x, v.z*x};
}

typedef struct {
	vec3 position;
	float v_w;
	float v_h;
	float distance;
} camera;


typedef struct {
	float r;
	float g;
	float b;
} color;

typedef struct {
	vec3 position;
	float radius;
	color color;
	float specular; // How shiny
	float reflective; // How reflective
} sphere;

void InitialiseImage(color image[IMG_H][IMG_W], color c) {
	for (int y = 0; y < IMG_H; ++y)
		for (int x = 0; x < IMG_W; ++x) {
			image[y][x] = c;
		}
}
color SetColor(unsigned char r, unsigned char g, unsigned char b) {
	// Set color as RGB Value
	color c = {r/255., g/255., b/255.};
	return c;
}
vec3 Vec3(float x, float y, float z) {
	vec3 v = {x,y,z};
	return v;
}
void PutPixel(int x, int y, color c, color image[IMG_H][IMG_W]) {
	// Convert x and y into screen coordinates
	int s_x = (IMG_W/2) + x;
	int s_y = (IMG_H/2) - y;
	if (s_x < 0 || s_x >= IMG_W) return;
	if (s_y < 0 || s_y >= IMG_H) return;
	image[s_y][s_x] = c;
}

vec3 CanvasToViewport(float x, float y, camera cam) {
	vec3 p = {x*cam.v_w/IMG_W, y*cam.v_h/IMG_H, cam.distance};
	return p;
}

sphere Sphere(vec3 position, float r, color col, float specular, float reflective) {
	sphere s = {position, r, col, specular, reflective};
	return s;
}
camera SetCamera(vec3 position, float v_w, float v_h, float distance) {
	camera c = {position, v_w, v_h, distance};
	return c;
}
void PrintVector(vec3 v) {
	printf("x=%f y=%f z=%f\n", v.x,v.y,v.z);
}

vec3 ReflectRay(vec3 R, vec3 N) {
	return Sub(MultScalar(N, 2*Dot(N,R)), R);
}

vec2 IntersectRaySphere(vec3 Origin, vec3 D, sphere sphere) {
	float r = sphere.radius;
	vec3 CO = Sub(Origin, sphere.position);
	float a = Dot(D, D);
	float b = 2*Dot(CO, D);
	float c = Dot(CO,CO) - (r*r);

	float discriminant = (b*b) - (4*a*c);
	// Intersections
	vec2 t = {T_MAX, T_MAX};
	if (discriminant < 0) {
		return t;
	} else {
		// t1
		t.x = (-b + sqrt(discriminant)) / (2*a);
		t.y = (-b - sqrt(discriminant)) / (2*a);
		return t;
	}
}

// Structure representing intersection between a sphere and its t value
struct SphereIntersection {
	unsigned char intersection;
	float t;
	sphere sphere;
};

struct SphereIntersection ClosestIntersection(vec3 O, vec3 D, float t_min, float t_max, sphere spheres[SPHERE_NUM]) {
	float closest_t = T_MAX;
	sphere* closest_sphere = NULL;
	for (int i = 0; i < SPHERE_NUM; ++i) {
		vec2 t = IntersectRaySphere(O, D, spheres[i]);
		float t1 = t.x;
		float t2 = t.y;
		if ((t1 > t_min && t1 < t_max) && (t1 < closest_t)) {
			closest_t = t1;
			closest_sphere = &spheres[i];
		}
		if ((t2 > t_min && t2 < t_max) && (t2 < closest_t)) {
			closest_t = t2;
			closest_sphere = &spheres[i];
		}
	}
	if (closest_sphere != NULL) {
		return (struct SphereIntersection){.intersection = 1, .t = closest_t, .sphere = *closest_sphere};
	} else {
		return (struct SphereIntersection){.intersection = 0};
	}
}


typedef enum {
	LIGHT_AMBIENT,
	LIGHT_POINT,
	LIGHT_DIRECTIONAL,
} LIGHT_T;

typedef struct {
	LIGHT_T type;
	float intensity;
	vec3 position;
	vec3 direction;
} light;

light DirectionalLight(float intensity, vec3 direction) {
	return (light){LIGHT_DIRECTIONAL, intensity, (vec3){0,0,0}, direction};
}

light PointLight(float intensity, vec3 position) {
	return (light){LIGHT_POINT, intensity, position, (vec3){0,0,0}};
}

light AmbientLight(float intensity) {
	return (light){LIGHT_AMBIENT, intensity, (vec3){0,0,0}, (vec3){0,0,0}};
}
color ColorMult(color c, float x) {
	return (color){((float)(c.r * x)), (float)(c.g * x), (float)(c.b * x)};
}
color ColorSum(color a, color b) {
	return (color){
		.r = (a.r + b.r),
		.g = (a.g + b.g),
		.b = (a.b + b.b),
	};
}

float clamp(float x, float lo, float hi) {
	if (x > hi) return hi;
	if (x < lo) return lo;
	return x;
}
float ComputeLighting(light lights[LIGHT_NUM], sphere spheres[SPHERE_NUM], vec3 D, vec3 Normal, vec3 V, float specular) {
	// Put scene lights here (hacky)
	//printf("%f\n", Mag(Normal));
	float mag_v = Mag(V);
	float intensity = 0;
	for (int i = 0; i < LIGHT_NUM; ++i) {
		light* light = &lights[i];
		vec3 L;
		if (light->type == LIGHT_AMBIENT) {
			intensity += light->intensity;
		} else {
			if (light->type == LIGHT_POINT) {
				L = Sub(light->position, D);
			} else {
				L = light->direction;
			}

			// Shadows
			struct SphereIntersection shadow_intersect = ClosestIntersection(D, L, 0.001, T_MAX, spheres);
			if (shadow_intersect.intersection == 1) {
				continue; // Shadow Magic
			}
			// Diffuse
			float light_dot_normal = Dot(Normal, L);
			if (light_dot_normal > 0.0) {
				//printf("%f\n", light->intensity * (light_dot_normal/(Mag(Normal)*Vec3Mag(L))));
				intensity += light->intensity * (light_dot_normal/(Mag(Normal)*Mag(L)));
			}

			// Specular
			if (specular != -1) {
				// Specularity involves the reflection of light and the difference of that angle between the viewpoint you get me
				// V is equal to Negative of D (pointing to the camera)
				vec3 Reflect = ReflectRay(L, Normal);
				float R_dot_V = Dot(Reflect, V);
				if (R_dot_V > 0) {
					float x = light->intensity * pow(R_dot_V/(Mag(Reflect)*mag_v), specular);
					if (x > 0) {
						intensity += light->intensity * pow(R_dot_V/(Mag(Reflect)*mag_v), specular);
					}
				}
			}
		}
	}
	return intensity;
}


color TraceRay(color default_bg, vec3 Origin, vec3 D, sphere spheres[SPHERE_NUM], light lights[LIGHT_NUM], float t_min, float t_max, int recursion_depth) {
	struct SphereIntersection intersect = ClosestIntersection(Origin, D, t_min, t_max, spheres);
	sphere *closest_sphere = NULL;
	if (intersect.intersection == 1) {
		closest_sphere = &intersect.sphere;
		float closest_t = intersect.t;
		vec3 intersection_point = Add(Origin,MultScalar(D,closest_t));
		vec3 Normal = Sub(intersection_point, closest_sphere->position);
		Normal = MultScalar(Normal, 1./Mag(Normal));
		// Get local color
		color localColor = ColorMult(closest_sphere->color, ComputeLighting(lights, spheres, intersection_point, Normal, MultScalar(D, -1.0), closest_sphere->specular));
		// If we hit recursion limit or the object is not reflective, we're done
		float reflect = closest_sphere->reflective;
		if (recursion_depth <= 0 || reflect <= 0) {
			return localColor;
		}
		vec3 R = ReflectRay(MultScalar(D,-1), Normal);
		color reflectedColor = TraceRay(default_bg, intersection_point, R, spheres, lights, 0.001, T_MAX, recursion_depth - 1);
		return ColorSum(ColorMult(localColor, (1 - reflect)), ColorMult(reflectedColor, reflect));

	}
	return default_bg;
}
int main() {
	color img[IMG_H][IMG_W];
	color default_bg = {0,0,0};
	InitialiseImage(img, default_bg);
	vec3 origin = {0,0,-5};
	camera Camera = SetCamera(origin, 1, 1, 1);
	light lights[LIGHT_NUM];
	lights[0] = AmbientLight(0.2);
	lights[1] = PointLight(0.6, Vec3(2,1,0));
	lights[2] = DirectionalLight(0.2, Vec3(1,4,4));
	sphere spheres[SPHERE_NUM];
	spheres[0] = Sphere(Vec3(0,-1,3), 1, SetColor(255,0,0), 500.0, 0.8);
	spheres[1] = Sphere(Vec3(2,0,4), 1, SetColor(0,0,255), 500.0, 0.3);
	spheres[2] = Sphere(Vec3(-2,0,4), 1, SetColor(0,255,0), 10.0, 0.4);
	spheres[3] = Sphere(Vec3(0,-5001,4), 5000, SetColor(255,255,0), 1000.0, 0);

	for (int x = -IMG_W/2; x < IMG_W/2; ++x) {
		for (int y = -IMG_H/2; y < IMG_H/2; ++y) {
			vec3 D = CanvasToViewport(x, y, Camera);
			//printf("%f %f %f\n", D.x, D.y, D.z);
			color C = TraceRay(default_bg, Camera.position, D, spheres, lights, T_MIN, T_MAX, RECURSION_LIMIT);
			PutPixel(x, y, C, img);
		}
	}
	// Output image data to file
	FILE *o_img;
	o_img = fopen(IMG_NAME, "wb");
	if (!o_img) {
		printf("Error opening file %s!\n", IMG_NAME);
		return 1;
	}

	/* Print image and close file */
	fprintf(o_img, "P6 %d %d 255 ", IMG_W, IMG_H); // Initialise PPM File output
	for (int y = 0; y < IMG_H; ++y) {
		for (int x = 0; x < IMG_W; ++x) {
			color c = img[y][x];
			unsigned char r = (unsigned char)(clamp(c.r,0,1)*255);
			unsigned char g = (unsigned char)(clamp(c.g,0,1)*255);
			unsigned char b = (unsigned char)(clamp(c.b,0,1)*255);
			fwrite(&r, 1, 1, o_img);
			fwrite(&g, 1, 1, o_img);
			fwrite(&b, 1, 1, o_img);
		}
	}
	fclose(o_img);
	return 0;
}
