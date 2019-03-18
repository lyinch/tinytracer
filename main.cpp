#include <vector>
#include <fstream>
#include "geometry.h"


struct Material {
    Material(const Vec3f &a, const Vec3f &color, const float &spec) : albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : albedo(1,0,0), diffuse_color(), specular_exponent() {}
    Vec3f diffuse_color;
    Vec3f albedo; //diffuse reflection parameter
    float specular_exponent;

};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Sphere{
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float &r,const Material &m): center(c), radius(r), material(m){};
#if 0
    // Geometric solution
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const{
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
#else
    // Analytic solution
    bool ray_intersect(const Vec3f &orig,  const Vec3f &dir, float &t0) const{
        Vec3f L = orig - center;
        float a = dir*dir;
        float b =  2*(dir*L);
        float c = L*L - radius*radius;
        //std::cout << L << "--- " << radius << std::endl;
        //std::cout << c << std::endl;
        // solve quadratic equation: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

        float discr = b*b-4*a*c;
        if (discr < 0) return false;
        float t1;
        if(discr == 0){
            t0 = t1 = -0.5f*b/a;
        }
        else{
            float q = (b > 0) ?
                      -0.5f * (b + sqrtf(discr)) :
                      -0.5f * (b - sqrtf(discr));
            t0 = q / a;
            t1 = c / q;
        }
        if (t0 > t1) std::swap(t0, t1);

        if (t0 < 0) {
            t0 = t1; // if t0 is negative, let's use t1 instead
            if (t0 < 0) return false; // both t0 and t1 are negative
        }
        return true;

    }
#endif
};

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for(auto sphere:spheres){
        float dist_i;
        if(sphere.ray_intersect(orig,dir,dist_i) && dist_i < spheres_dist){
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit-sphere.center).normalize();
            material = sphere.material;
        }
    }
    return spheres_dist<1000;

}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir,const std::vector<Sphere> &spheres, const std::vector<Light> &lights,size_t depth=0) {
    Vec3f point, N;
    Material material;
    float diffuse_light_intensity = 0,specular_light_intensity = 0;;

    if(depth>4 || !scene_intersect(orig,dir,spheres,point,N, material)){
        return {0.2, 0.7, 0.8}; // background color
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);

    for(auto light:lights){
        Vec3f light_dir = (light.position-point).normalize();

        float light_distance = (light.position - point).norm();
        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // create the origin for the shadow ray. Move the ray a bit to not intersect with the current sphere
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;

        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += light.intensity*std::max(0.f,light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*light.intensity;
    }

    return material.diffuse_color*diffuse_light_intensity*material.albedo[0] + Vec3f(1.f, 1.f, 1.f)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2];
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights){
    const int width = 1024;
    const int height = 768;
    const auto fov      = static_cast<const float>(M_PI / 2.f);

    std::vector<Vec3f> framebuffer(width*height);
    for(size_t j = 0; j < height; j++){
        for(size_t i = 0; i < width; i++){
            float x = static_cast<float>((2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height);
            float y = static_cast<float>(-(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.));
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, spheres, lights);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < height * width; ++i) {
        //normalize value
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (int j = 0; j < 3; ++j) {
            ofs << (char)(255*std::max(0.f,std::min(1.f,framebuffer[i][j])));
        }
    }
}


int main() {
    Material      ivory(Vec3f(0.6,  0.3, 0.1), Vec3f(0.4, 0.4, 0.3),   50.);
    Material red_rubber(Vec3f(0.9,  0.1, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(Vec3f(0.0, 10.0, 0.8), Vec3f(1.0, 1.0, 1.0), 800.);

    std::vector<Sphere> spheres;
    spheres.emplace_back(Vec3f(-3,    0,   -16), 2,ivory);
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,      mirror));
    spheres.emplace_back(Vec3f( 1.5, -0.f, -18), 3,red_rubber);
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,     mirror));


    std::vector<Light>  lights;
    lights.emplace_back(Vec3f(-20, 20,  20), 1.5);
    lights.emplace_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.emplace_back(Light(Vec3f( 30, 20,  30), 1.7));

    render(spheres, lights);



    return 0;
}