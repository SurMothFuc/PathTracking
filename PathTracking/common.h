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
T Clamp(T x, T min, T max)//�ضϺ���
{
    if (x > max)
        return max;
    if (x < min)
        return min;
    return x;
}
class Material //������
{
public:
    bool isEmissive = false;        // �Ƿ񷢹�
    vec3 normal = vec3(0, 0, 0);    // ������
    vec3 color = vec3(0, 0, 0);     // ��ɫ
    double specularRate = 0.0f;      // �����ռ��
};
class HitResult//�󽻵ĸ����� ���潻��Ľ��
{
public:
    bool isHit = false;             // �Ƿ�����
    double distance = 0.0f;         // �뽻��ľ���
    vec3 hitPoint = vec3(0, 0, 0);  // �������е�
    Material material;              // ���е�ı������
};
class Ray
{
    public:
    vec3 startPoint = vec3(0, 0, 0);    // ���
    vec3 direction = vec3(0, 0, 0);     // ����
};
class Shape
{
public:
    Shape() {}
    virtual HitResult intersect(Ray ray) { return HitResult(); }//���м̳л���ĺ���Ӧ��ʵ���󽻵ĺ��� �������󽻵Ľ��
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
    vec3 p1, p2, p3;    // ������
    Material material;  // ����

    // �������
    HitResult intersect(Ray ray)
    {
        HitResult res;
        vec3 S = ray.startPoint;
        vec3 d = ray.direction;
        vec3 N = material.normal;

        if (N.dot(d) > 0.0f)//˵�������Ǵ������α���������� ��Ҫ��ת������
            N = -N;

        if (abs(N.dot(d)) < 0.00001f)//˵����������������ƽ�е� ��ôֱ�ӷ���
            return res;

        float t = (N.dot(p1) - S.dot(N)) / d.dot(N);
        if (t < 0.0005f)//��������������� 0.0005�Ƿ�ֹ��������
            return res;
        vec3 P = S + d * t;

        vec3 c1 = (p2 - p1).cross(P - p1);
        vec3 c2 = (p3 - p2).cross(P - p2);
        vec3 c3 = (p1 - p3).cross( P - p3);
        vec3 n = material.normal;
        if (c1.dot(n) < 0 || c2.dot(n) < 0 || c3.dot(n) < 0)//�㲻���������ڲ������
            return res;

        res.isHit = true;
        res.distance = t;
        res.hitPoint = P;
        res.material = material;
        res.material.normal = N;

        return res;
    };
};
uniform_real_distribution<> dis(0.0, 1.0); //�Ӿ��ȷֲ�����������ĸ�����
random_device rd;//����һ��Ϊ�����
mt19937  rand_num(rd());//��rd��Ϊ����
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
vec3 randomDirection(vec3 n)//�����ڰ����ϵ��������
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