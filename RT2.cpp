#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

int Height = 500;
int Width = 500;

struct Vec3 {
	double x;
	double y;
	double z;
	Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	Vec3 operator +(Vec3 v) { return Vec3(x + v.x, y + v.y, z + v.z); }
	Vec3 operator -(Vec3 v) { return Vec3(x - v.x, y - v.y, z - v.z); }
	Vec3 operator *(double d) { return Vec3(x*d, y*d, z*d); }
	Vec3 operator /(double d) { return Vec3(x / d, y / d, z / d); }
	Vec3 normalize() {
		double mg = sqrt(x*x + y*y + z*z);
		return Vec3(x / mg, y / mg, z / mg);
	}
};

Vec3 white(255, 255, 255);
Vec3 black(0, 0, 0);
Vec3 blue(0, 0, 255);
Vec3 red(255, 0, 0);
Vec3 green(0, 255, 0);
Vec3 yellow(255, 255, 0);
Vec3 cyan(0, 255, 255);

Vec3 pyr1(0, 0, 0);
Vec3 pyr2(0, 0, 0);
Vec3 pyr3(0, 0, 0);
Vec3 pyr4(0, 0, 0);

inline double dot(Vec3 a, Vec3 b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline Vec3 cross(Vec3 a, Vec3 b) {
	return (Vec3((a.y*b.z - b.y*a.z), (b.x*a.z - a.x*b.z), (a.x*b.y - b.x*a.y)));
}

Vec3 Make_Pyramid1(Vec3 c, double l, double rot)
{
	pyr1 = c;
	return pyr1;
}

Vec3 Make_Pyramid2(Vec3 c, double l, double rot)
{
	pyr2 = c + Vec3((l*cos(rot)), 0, (-l*sin(rot)));
	return pyr2;
}

Vec3 Make_Pyramid3(Vec3 c, double l, double rot)
{
	pyr3 = c + Vec3((((l / 2)*cos(rot)) - ((sqrt(3.0)*l / 2)*sin(rot))), 0, ((-(sqrt(3.0)*l / 2)*cos(rot)) - (l*sin(rot) / 2)));
	return pyr3;
}

Vec3 Make_Pyramid4(Vec3 c, double l, double rot)
{
	pyr4 = c + Vec3((((l / 2)*cos(rot)) - ((l / (2 * sqrt(3.0)))*sin(rot))), (l*sqrt(2.0) / sqrt(3.0)), (((-l / (2 * sqrt(3.0)))*cos(rot)) - (l*sin(rot) / 2)));
	return pyr4;
}

struct Ray {
	Vec3 o, d;
	Ray(Vec3 o, Vec3 d) : o(o), d(d) {}
};

struct Sphere {
	Vec3 c;
	double r, t_v;
	Sphere(Vec3 c, double r) : c(c), r(r) {}
	Vec3 getNormal(Vec3 pi) { return (pi - c) / r; }
	bool intersect(Ray ray, double t) {
		Vec3 o = ray.o;
		Vec3 d = ray.d;
		Vec3 oc = o - c;
		double b = 2 * dot(oc, d);
		double c = dot(oc, oc) - r*r;
		double disc = b*b - 4 * c;
		if (disc < 0) return false;
		disc = sqrt(disc);
		double t0 = (-b - disc) / 2;
		double t1 = (-b + disc) / 2;
		if (t0 > 0.001 && t1 > 0.001)
			t = (t0 < t1) ? t0 : t1;
		else if (t1 > 0.001)
			t = t1;
		else if (t0 > 0.001)
			t = t0;
		else
			t = -INT_MAX;
		t_v = t;
		if (t_v >= 0)
			return true;
		else
			return false;
	}
	double t_val() {
		if (t_v >= 0)
			return t_v;
		else
			return INT_MAX;
	}
};

struct Rectangle {
	Vec3 p1, p2, p3, p4;
	double t_v;
	Rectangle(Vec3 p1, Vec3 p2, Vec3 p3, Vec3 p4) : p1(p1), p2(p2), p3(p3), p4(p4) {}
	Vec3 getNormal() { return cross((p2 - p1), (p4 - p1)).normalize(); }
	bool intersect(Ray ray, double t) {
		Vec3 o = ray.o;
		Vec3 d = ray.d;
		Vec3 n = getNormal();
		Vec3 c = ((p1 + p2 + p3 + p4) / 4);
		t = (dot((p1 - o), n)) / (dot(d, n));
		Vec3 p = o + d*t;
		t_v = t;
		if (t_v>0 && ((dot((p - ((p1 + p2) / 2)), (c - ((p1 + p2) / 2))) > 0) && (dot((p - ((p2 + p3) / 2)), (c - ((p2 + p3) / 2))) > 0) && (dot((p - ((p3 + p4) / 2)), (c - ((p3 + p4) / 2))) > 0) && (dot((p - ((p1 + p4) / 2)), (c - ((p1 + p4) / 2))) > 0)))
			return true;
		else
			return false;
	}
	double t_val() {
		if (t_v > 0)
			return t_v;
		else
			return INT_MAX;
	}
};

struct Pyramid {
	Vec3 p1, p2, p3;
	double t_v;
	Pyramid(Vec3 p1, Vec3 p2, Vec3 p3) : p1(p1), p2(p2), p3(p3) {}
	Vec3 getNormal() { return cross((p2 - p1), (p3 - p1)).normalize(); }
	bool intersect(Ray ray, double t) {
		Vec3 o = ray.o;
		Vec3 d = ray.d;
		Vec3 n = getNormal();
		Vec3 c = ((p1 + p2 + p3) / 3);
		t = (dot((p1 - o), n)) / (dot(d, n));
		Vec3 p = o + d*t;
		t_v = t;
		if (t_v > 0.0001 && ((dot((p - ((p1 + p2) / 2)), (c - ((p1 + p2) / 2))) > 0) && (dot((p - ((p2 + p3) / 2)), (c - ((p2 + p3) / 2))) > 0) && (dot((p - ((p3 + p1) / 2)), (c - ((p3 + p1) / 2))) > 0)))
		{
			return true;
		}
		else
			return false;
	}
	double t_val() {
		if (t_v > 0.0001)
			return t_v;
		else
			return INT_MAX;
	}
};

Sphere sphere(Vec3(275, 90, -150), 90);
Sphere light(Vec3(Width*0.5, Height - 50, -250), 1);
Rectangle ground(Vec3(0, 0, 0), Vec3(500, 0, 0), Vec3(500, 0, -500), Vec3(0, 0, -500));
Rectangle wall1(Vec3(0, 0, 0), Vec3(0, 0, -500), Vec3(0, 500, -500), Vec3(0, 500, 0));
Rectangle wall2(Vec3(500, 0, 0), Vec3(500, 500, 0), Vec3(500, 500, -500), Vec3(500, 0, -500));
Rectangle wall3(Vec3(0, 0, -500), Vec3(500, 0, -500), Vec3(500, 500, -500), Vec3(0, 500, -500));
Rectangle ceiling(Vec3(0, 500, 0), Vec3(0, 500, -500), Vec3(500, 500, -500), Vec3(500, 500, 0));


Vec3 pyr_origin(100, 0, -25);
double pyr_len = 170;
double pyr_rotation = 3.1415 / 3;
Pyramid pyr_wall1(Make_Pyramid1(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid2(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid4(pyr_origin, pyr_len, pyr_rotation));
Pyramid pyr_wall2(Make_Pyramid2(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid3(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid4(pyr_origin, pyr_len, pyr_rotation));
Pyramid pyr_wall3(Make_Pyramid3(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid1(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid4(pyr_origin, pyr_len, pyr_rotation));
Pyramid pyr_wall4(Make_Pyramid1(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid3(pyr_origin, pyr_len, pyr_rotation), Make_Pyramid2(pyr_origin, pyr_len, pyr_rotation));

Vec3 c_origin(120, 0, -260);
double c_len = 150;
double c_rotation = 0;
Rectangle c_wall1(c_origin + Vec3(0, 0, 0), c_origin + Vec3(c_len, 0, 0), c_origin + Vec3(c_len, c_len, 0), c_origin + Vec3(0, c_len, 0));
Rectangle c_wall2(c_origin + Vec3(c_len, 0, 0), c_origin + Vec3(c_len, 0, -c_len), c_origin + Vec3(c_len, c_len, -c_len), c_origin + Vec3(c_len, c_len, 0));
Rectangle c_wall3(c_origin + Vec3(0, 0, 0), c_origin + Vec3(0, c_len, 0), c_origin + Vec3(0, c_len, -c_len), c_origin + Vec3(0, 0, -c_len));
Rectangle c_wall4(c_origin + Vec3(c_len, 0, -c_len), c_origin + Vec3(0, 0, -c_len), c_origin + Vec3(0, c_len, -c_len), c_origin + Vec3(c_len, c_len, -c_len));
Rectangle c_ceiling(c_origin + Vec3(0, c_len, 0), c_origin + Vec3(c_len, c_len, 0), c_origin + Vec3(c_len, c_len, -c_len), c_origin + Vec3(0, c_len, -c_len));

Vec3 rayTrace(Ray ray)
{
	double arr[15];
	double t;
	int flag = 0;
	for (int i = 0; i < 15; i++)
	{
		arr[i] = INT_MAX;
	}

	if (sphere.intersect(ray, t)) {
		arr[0] = sphere.t_val();
		flag = 1;
	}
	if (ground.intersect(ray, t)) {
		arr[1] = ground.t_val();
		flag = 1;
	}
	if (wall1.intersect(ray, t)) {
		arr[2] = wall1.t_val();
		flag = 1;
	}
	if (wall2.intersect(ray, t)) {
		arr[3] = wall2.t_val();
		flag = 1;
	}
	if (wall3.intersect(ray, t)) {
		arr[4] = wall3.t_val();
		flag = 1;
	}
	if (ceiling.intersect(ray, t)) {
		arr[5] = ceiling.t_val();
		flag = 1;
	}
	if (pyr_wall1.intersect(ray, t)) {
		arr[6] = pyr_wall1.t_val();
		flag = 1;
	}
	if (pyr_wall2.intersect(ray, t)) {
		arr[7] = pyr_wall2.t_val();
		flag = 1;
	}
	if (pyr_wall3.intersect(ray, t)) {
		arr[8] = pyr_wall3.t_val();
		flag = 1;
	}
	if (pyr_wall4.intersect(ray, t)) {
		arr[9] = pyr_wall4.t_val();
		flag = 1;
	}
	if (c_wall1.intersect(ray, t)) {
		arr[10] = c_wall1.t_val();
		flag = 1;
	}
	if (c_wall2.intersect(ray, t)) {
		arr[11] = c_wall2.t_val();
		flag = 1;
	}
	if (c_wall3.intersect(ray, t)) {
		arr[12] = c_wall3.t_val();
		flag = 1;
	}
	if (c_wall4.intersect(ray, t)) {
		arr[13] = c_wall4.t_val();
		flag = 1;
	}
	if (c_ceiling.intersect(ray, t)) {
		arr[14] = c_ceiling.t_val();
		flag = 1;
	}

	if (flag == 0)
		return black;

	int min_ind = 0;

	for (int i = 0; i < 15; i++)
	{
		if (arr[min_ind] - arr[i] > double(0.0))
			min_ind = i;
	}

	if (min_ind == 0)
	{
		Vec3 pi = ray.o + (ray.d * arr[0]);
		Vec3 L = light.c - pi;
		Vec3 N = sphere.getNormal(pi);
		double dt = dot(L.normalize(), N.normalize());
		Vec3 Per = (ray.d - N * dot(ray.d, N)).normalize();
		Vec3 Nor = (Vec3(0.0, 0.0, 0.0) - N).normalize();
		double decide = dot(ray.d, Nor);
		double cosi = abs(decide);
		double sini = sqrt(1 - cosi*cosi);
		if (sini == double(0.0))
			sini = 0.000001;
		if (cosi == double(0.0))
			cosi = 0.000001;
		Ray ray3(pi, ray.d - N * 2 * dot(ray.d, N));
		if (decide > 0)
		{
			double sinr = (4.0 / 5.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			Vec3 refract(Nor*cosr+Per*sinr);
			Ray ray2(pi,refract);
			return rayTrace(ray2)*0.9 + rayTrace(ray3)*0.1;
		}
		else
		{
			double sinr = (5.0 / 4.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			Vec3 refract(Nor*-cosr + Per*sinr);
			Ray ray2(pi, refract);
			return rayTrace(ray2);
		}
	}
	else if (min_ind == 1)
	{
		Vec3 pi = ray.o + (ray.d * arr[1]);
		Vec3 L = light.c - pi;
		Vec3 N = ground.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		Ray ray2(pi, ray.d - N * 2 * dot(ray.d, N));
		if (c_wall1.intersect(ray_dummy, t_dummy) || c_wall2.intersect(ray_dummy, t_dummy) || c_wall3.intersect(ray_dummy, t_dummy) || c_wall4.intersect(ray_dummy, t_dummy) || c_ceiling.intersect(ray_dummy, t_dummy))
		{
			if ((int(pi.x / 50)) % 2 != (int(-pi.z / 50)) % 2)
				return (rayTrace(ray2)*0.25 + white*(0.5 + dt*0.5)*0.75)*0.2;
			else
				return (rayTrace(ray2)*0.25 + black*(0.5 + dt*0.5)*0.75)*0.2;
		}
		else
		{
			if ((int(pi.x / 50)) % 2 != (int(-pi.z / 50)) % 2)
				return rayTrace(ray2)*0.25 + white*(0.5 + dt*0.5)*0.75;
			else
				return rayTrace(ray2)*0.25 + black*(0.5 + dt*0.5)*0.75;
		}
	}
	else if (min_ind == 2)
	{
		Vec3 pi = ray.o + (ray.d * arr[2]);
		Vec3 L = light.c - pi;
		Vec3 N = wall1.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy) || c_wall1.intersect(ray_dummy, t_dummy) || c_wall2.intersect(ray_dummy, t_dummy) || c_wall3.intersect(ray_dummy, t_dummy) || c_wall4.intersect(ray_dummy, t_dummy) || c_ceiling.intersect(ray_dummy, t_dummy))
		{
			return blue*(dt)*0.2;
		}
		else
		{
			return blue*(dt);
		}
	}
	else if (min_ind == 3)
	{
		Vec3 pi = ray.o + (ray.d * arr[3]);
		Vec3 L = light.c - pi;
		Vec3 N = wall2.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy) || c_wall1.intersect(ray_dummy, t_dummy) || c_wall2.intersect(ray_dummy, t_dummy) || c_wall3.intersect(ray_dummy, t_dummy) || c_wall4.intersect(ray_dummy, t_dummy) || c_ceiling.intersect(ray_dummy, t_dummy))
		{
			return green*(dt)*0.2;
		}
		else
		{
			return green*(dt);
		}
	}
	else if (min_ind == 4)
	{
		Vec3 pi = ray.o + (ray.d * arr[4]);
		Vec3 L = light.c - pi;
		Vec3 N = wall3.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy) || c_wall1.intersect(ray_dummy, t_dummy) || c_wall2.intersect(ray_dummy, t_dummy) || c_wall3.intersect(ray_dummy, t_dummy) || c_wall4.intersect(ray_dummy, t_dummy) || c_ceiling.intersect(ray_dummy, t_dummy))
		{
			return yellow*(0.8 + 0.2*dt)*0.2;
		}
		else
		{
			return yellow*0.5*(0.8 + 0.2*dt);
		}
	}
	else if (min_ind == 5)
	{
		Vec3 pi = ray.o + (ray.d * arr[5]);
		Vec3 L = light.c - pi;
		Vec3 N = ceiling.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		if (pi.x >= 150 && pi.x <= 350 && pi.z <= -150 && pi.z >= -350)
			return white;
		return white*0.5*dt;
	}
	else if (min_ind == 6)
	{
		Vec3 pi = ray.o + (ray.d * arr[6]);
		Vec3 L = light.c - pi;
		Vec3 N = pyr_wall1.getNormal();
		double dt = dot(L.normalize(), N.normalize());
		Vec3 Per = (ray.d - N * dot(ray.d, N)).normalize();
		Vec3 Nor = (Vec3(0.0, 0.0, 0.0) - N).normalize();
		double decide = dot(ray.d, Nor);
		double cosi = abs(decide);
		double sini = sqrt(1 - cosi*cosi);
		if (sini == double(0.0))
			sini = 0.000001;
		if (cosi == double(0.0))
			cosi = 0.000001;
		Ray ray3(pi, ray.d - N * 2 * dot(ray.d, N));
		if (decide > 0)
		{
			double sinr = (4.0 / 5.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			Vec3 refract(Nor*cosr + Per*sinr);
			Ray ray2(pi, refract);
			return rayTrace(ray2)*0.9 + rayTrace(ray3)*0.1;
		}
		else
		{
			double sinr = (5.0 / 4.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			if (sini < (4.0 / 5.0))
			{
				Vec3 refract(Nor*-cosr + Per*sinr);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
			else
			{
				Vec3 refract(Nor*cosi + Per*sini);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
		}
	}
	else if (min_ind == 7)
	{
		Vec3 pi = ray.o + (ray.d * arr[7]);
		Vec3 L = light.c - pi;
		Vec3 N = pyr_wall2.getNormal();
		double dt = dot(L.normalize(), N.normalize());
		Vec3 Per = (ray.d - N * dot(ray.d, N)).normalize();
		Vec3 Nor = (Vec3(0.0, 0.0, 0.0) - N).normalize();
		double decide = dot(ray.d, Nor);
		double cosi = abs(decide);
		double sini = sqrt(1 - cosi*cosi);
		if (sini == double(0.0))
			sini = 0.000001;
		if (cosi == double(0.0))
			cosi = 0.000001;
		Ray ray3(pi, ray.d - N * 2 * dot(ray.d, N));
		if (decide > 0)
		{
			double sinr = (4.0 / 5.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			Vec3 refract(Nor*cosr + Per*sinr);
			Ray ray2(pi, refract);
			return rayTrace(ray2)*0.9 + rayTrace(ray3)*0.1;
		}
		else
		{
			double sinr = (5.0 / 4.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			if (sini < (4.0 / 5.0))
			{
				Vec3 refract(Nor*-cosr + Per*sinr);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
			else
			{
				Vec3 refract(Nor*cosi + Per*sini);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
		}
	}
	else if (min_ind == 8)
	{
		Vec3 pi = ray.o + (ray.d * arr[8]);
		Vec3 L = light.c - pi;
		Vec3 N = pyr_wall3.getNormal();
		double dt = dot(L.normalize(), N.normalize());
		Vec3 Per = (ray.d - N * dot(ray.d, N)).normalize();
		Vec3 Nor = (Vec3(0.0, 0.0, 0.0) - N).normalize();
		double decide = dot(ray.d, Nor);
		double cosi = abs(decide);
		double sini = sqrt(1 - cosi*cosi);
		if (sini == double(0.0))
			sini = 0.000001;
		if (cosi == double(0.0))
			cosi = 0.000001;
		Ray ray3(pi, ray.d - N * 2 * dot(ray.d, N));
		if (decide > 0)
		{
			double sinr = (4.0 / 5.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			Vec3 refract(Nor*cosr + Per*sinr);
			Ray ray2(pi, refract);
			return rayTrace(ray2)*0.9 + rayTrace(ray3)*0.1;
		}
		else
		{
			double sinr = (5.0 / 4.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			if (sini < (4.0 / 5.0))
			{
				Vec3 refract(Nor*-cosr + Per*sinr);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
			else
			{
				Vec3 refract(Nor*cosi + Per*sini);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
		}
	}
	else if (min_ind == 9)
	{
		Vec3 pi = ray.o + (ray.d * arr[9]);
		Vec3 L = light.c - pi;
		Vec3 N = pyr_wall4.getNormal();
		double dt = dot(L.normalize(), N.normalize());
		Vec3 Per = (ray.d - N * dot(ray.d, N)).normalize();
		Vec3 Nor = (Vec3(0.0, 0.0, 0.0) - N).normalize();
		double decide = dot(ray.d, Nor);
		double cosi = abs(decide);
		double sini = sqrt(1 - cosi*cosi);
		if (sini == double(0.0))
			sini = 0.000001;
		if (cosi == double(0.0))
			cosi = 0.000001;
		Ray ray3(pi, ray.d - N * 2 * dot(ray.d, N));
		if (decide > 0)
		{
			double sinr = (4.0 / 5.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			Vec3 refract(Nor*cosr + Per*sinr);
			Ray ray2(pi, refract);
			return rayTrace(ray2)*0.9 + rayTrace(ray3)*0.1;
		}
		else
		{
			double sinr = (5.0 / 4.0)*sini;
			double cosr = sqrt(1 - sinr*sinr);
			if (sini < (4.0 / 5.0))
			{
				Vec3 refract(Nor*-cosr + Per*sinr);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
			else
			{
				Vec3 refract(Nor*cosi + Per*sini);
				Ray ray2(pi, refract);
				return rayTrace(ray2);
			}
		}
	}
	else if (min_ind == 10)
	{
		Vec3 pi = ray.o + (ray.d * arr[9]);
		Vec3 L = light.c - pi;
		Vec3 N = c_wall1.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		Ray ray2(pi, ray.d - N * 2 * dot(ray.d, N));
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy))
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7)*0.2;
		}
		else
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7);
		}
	}
	else if (min_ind == 11)
	{
		Vec3 pi = ray.o + (ray.d * arr[10]);
		Vec3 L = light.c - pi;
		Vec3 N = c_wall2.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		Ray ray2(pi, ray.d - N * 2 * dot(ray.d, N));
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy))
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7)*0.2;
		}
		else
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7);
		}
	}
	else if (min_ind == 12)
	{
		Vec3 pi = ray.o + (ray.d * arr[11]);
		Vec3 L = light.c - pi;
		Vec3 N = c_wall3.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		Ray ray2(pi, ray.d - N * 2 * dot(ray.d, N));
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy))
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7)*0.2;
		}
		else
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7);
		}
	}
	else if (min_ind == 13)
	{
		Vec3 pi = ray.o + (ray.d * arr[12]);
		Vec3 L = light.c - pi;
		Vec3 N = c_wall4.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		Ray ray2(pi, ray.d - N * 2 * dot(ray.d, N));
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy))
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7)*0.2;
		}
		else
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7);
		}
	}
	else if (min_ind == 14)
	{
		Vec3 pi = ray.o + (ray.d * arr[2]);
		Vec3 L = light.c - pi;
		Vec3 N = c_ceiling.getNormal();
		double dt = dot(L.normalize(), N.normalize()), t_dummy = 0;
		Ray ray_dummy(light.c, (pi - light.c).normalize());
		Ray ray2(pi, ray.d - N * 2 * dot(ray.d, N));
		if (sphere.intersect(ray_dummy, t_dummy) || pyr_wall1.intersect(ray_dummy, t_dummy) || pyr_wall2.intersect(ray_dummy, t_dummy) || pyr_wall3.intersect(ray_dummy, t_dummy))
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7)*0.2;
		}
		else
		{
			return (rayTrace(ray2)*0.3 + cyan*(0.5 + 0.5*dt)*0.7);
		}
	}
	else
		return black;
}

int main() {

	std::ofstream out("output2.ppm");
	out << "P3\n" << Width << ' ' << Height << ' ' << "255\n";
	Vec3 ve = cross(Vec3(1, 0, 0), Vec3(0, 1, 0));

	Vec3 pix_col(black);

	for (int y = Height - 1; y >= 0; --y) {
		for (int x = 0; x < Width; ++x) {
			Ray ray(Vec3(x, y, 0), (Vec3(x, y, 0) - Vec3(Width/2, Height/2, 400)).normalize());
			pix_col = rayTrace(ray);
			out << (int)pix_col.x << ' '
				<< (int)pix_col.y << ' '
				<< (int)pix_col.z << '\n';
		}
	}
}