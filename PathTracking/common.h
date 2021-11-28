#pragma once

#include <Eigen/Dense>
#include <svpng.inc>
#include <vector>
#include <omp.h>
#include <fstream>
#include <random>

#define vec3 Vector3f

using namespace std;
using namespace Eigen;

template<class T>
T Clamp(T x, T min, T max)//截断函数
{
    if (x > max)
        return max;
    if (x < min)
        return min;
    return x;
}
class Material //材质类
{
public:
    bool isEmissive = false;        // 是否发光
    vec3 normal = vec3(0, 0, 0);    // 法向量
    vec3 color = vec3(0, 0, 0);     // 颜色
    double specularRate = 0.0f;      // 反射光占比
};
class HitResult//求交的辅助类 保存交点的结果
{
public:
    bool isHit = false;             // 是否命中
    double distance = 0.0f;         // 与交点的距离
    vec3 hitPoint = vec3(0, 0, 0);  // 光线命中点
    Material material;              // 命中点的表面材质
};
class Ray
{
    public:
    vec3 startPoint = vec3(0, 0, 0);    // 起点
    vec3 direction = vec3(0, 0, 0);     // 方向
};
class Shape
{
public:
    Shape() {}
    virtual HitResult intersect(Ray ray) { return HitResult(); }//所有继承基类的函数应该实现求交的函数 并返回求交的结果
};
class Triangle : public Shape
{
public:
    Triangle() {}
    Triangle(Vector3f P1, Vector3f P2, Vector3f P3, Vector3f C)
    {
        p1 = P1, p2 = P2, p3 = P3;
        material.normal = ((p2 - p1).cross( p3 - p1)).normalized(); 
        material.color = C;
    }
    vec3 p1, p2, p3;    // 三顶点
    Material material;  // 材质

    // 与光线求交
    HitResult intersect(Ray ray)
    {
        HitResult res;
        vec3 S = ray.startPoint;
        vec3 d = ray.direction;
        vec3 N = material.normal;

        if (N.dot(d) > 0.0f)//说明光线是从三角形背面射过来的 需要翻转法向量
            N = -N;

        if (abs(N.dot(d)) < 0.00001f)//说明光线是与三角形平行的 那么直接返回
            return res;

        float t = (N.dot(p1) - S.dot(N)) / d.dot(N);
        if (t < 0.0005f)//三角形在摄像机后 0.0005是防止精度问题
            return res;
        vec3 P = S + d * t;

        vec3 c1 = (p2 - p1).cross(P - p1);
        vec3 c2 = (p3 - p2).cross(P - p2);
        vec3 c3 = (p1 - p3).cross( P - p3);
        vec3 n = material.normal;
        if (c1.dot(n) < 0 || c2.dot(n) < 0 || c3.dot(n) < 0)//点不在三角形内部的情况
            return res;

        res.isHit = true;
        res.distance = t;
        res.hitPoint = P;
        res.material = material;
        res.material.normal = N;

        return res;
    };
};
uniform_real_distribution<> dis(0.0, 1.0); //从均匀分布中生成随机的浮点数
random_device rd;//生成一个为随机数
mt19937  rand_num(rd());//把rd设为种子
double randf()
{
    return dis(rand_num);
}
vec3 randomVec3()
{
    vec3 d;
    do
    {
        d = 2.0f * vec3(randf(), randf(), randf()) - vec3(1, 1, 1);
    } while (d.dot(d) > 1.0);
    return d.normalized();
}
vec3 randomDirection(vec3 n)//生成在半球上的随机变量
{
    return (randomVec3() + n).normalized();
}

class ProcessBar {
public:
    int rate;
    int length;
    int sum;
    clock_t t;
    clock_t startt;
    std::string info;
    ProcessBar(int s, std::string information = "") :rate(0), length(50), info(information), sum(s) {}
    void setInfo(std::string sinfo) {
        info = sinfo;
    }
    void start() {
        startt = clock();
        t = 0;
        update(rate, true);
    }
    void end() {
        rate = sum;
        update(rate, true);
        printf("\n");
    }
    void update(int schedule, bool f = false) {
        int oldt = t;
        t = clock() - startt;
        if (1.0 * (schedule - rate) / sum < 0.001 && t - oldt < 1000 && !f) {
            t = oldt;
            return;
        }
        int lc, rc;
        rate = schedule;
        lc = (int)(1.0 * rate / sum * length);
        rc = length - lc;
        printf("\r%s", info.c_str());
        printf("[");
        for (int i = 0; i < lc; i++)
            printf(">");
        for (int i = 0; i < rc; i++)
            printf("=");
        printf("]");
        printf("%.1f%% ", 100.0 * rate / sum);
        int second = (int)(1.0 * t / 1000);
        if (second / 3600 > 99)
            printf("runtime:99:99:99 ");
        else
            printf("runtime:%02d:%02d:%02d ", second / 3600, second % 3600 / 60, second % 60);

        float avespeed = 1.0 * rate / t * 1000;
        int esecond = (int)((sum - rate) / avespeed);

        if (esecond / 3600 > 99 || t == 0)
            printf("eta:99:99:99");
        else
            printf("eta:%02d:%02d:%02d ", esecond / 3600, esecond % 3600 / 60, esecond % 60);
        printf("\r");
    }
};
vec3 reflect(vec3 I, vec3 N)
{
    return I - N * 2 * (I.transpose() * N);
}