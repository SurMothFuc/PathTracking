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
    double refractRate = 0.0f;      // 折射光占比
    double refractAngle = 1.0f;     // 折射率
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

    vec3 center;       // 中心 用于bvh构建
    Shape() {}
    virtual HitResult intersect(Ray ray) { return HitResult(); }//所有继承基类的函数应该实现求交的函数 并返回求交的结果
    virtual vec3 minAA() { return vec3(0.0, 0.0, 0.0); }
    virtual vec3 maxBB() { return vec3(0.0, 0.0, 0.0); }
};
class Sphere : public Shape
{
public:
    Sphere() {}
    Sphere(vec3 o, double r, vec3 c) { O = o; R = r; material.color = c; center = O; }
    vec3 O;             // 圆心
    double R;           // 半径
    Material material;  // 材质

    // 与光线求交
    HitResult intersect(Ray ray)
    {
        HitResult res;

        vec3 S = ray.startPoint;        // 射线起点
        vec3 d = ray.direction;         // 射线方向

        float OS = (O - S).norm();
        float SH = (O - S).dot(d);
        float OH = sqrt(pow(OS, 2) - pow(SH, 2));

        if (OH > R) return res; // OH大于半径则不相交

        float PH = sqrt(pow(R, 2) - pow(OH, 2));

        float t1 = SH - PH;
        float t2 = SH + PH;
        float t = (t1 < 0) ? (t2) : (t1);   // 最近距离
        vec3 P = S + t * d;     // 交点

        // 防止自己交自己
        if (fabs(t1) < 0.0005f || fabs(t2) < 0.0005f) return res;

        // 装填返回结果
        res.isHit = true;
        res.distance = t;
        res.hitPoint = P;
        res.material = material;
        res.material.normal = (P - O).normalized(); // 要返回正确的法向
        return res;
    }
    vec3 minAA() {
        return vec3(O.x()-R,O.y()-R,O.z()-R);
    }
    vec3 maxBB() {
        return vec3(O.x() + R, O.y() + R, O.z() + R);
    }
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
        center = (P1 + P2 + P3) / 3;
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
    vec3 minAA() { 
        return vec3(
            min(p1.x(), min(p2.x(), p3.x())),
            min(p1.y(), min(p2.y(), p3.y())),
            min(p1.z(), min(p2.z(), p3.z())));
    }
    vec3 maxBB() {
        return vec3(
            max(p1.x(), max(p2.x(), p3.x())),
            max(p1.y(), max(p2.y(), p3.y())),
            max(p1.z(), max(p2.z(), p3.z())));
    }




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
bool refract(vec3 I, vec3 N, float eta,vec3& refracted)
{
    float a = I.transpose() * N;
    float cosi = Clamp(a, -1.0f, 1.0f);
    vec3 n = N;
    if (cosi < 0) { cosi = -cosi; }
    else {eta=1.0/eta; n = -N; }
    float k = 1 - eta * eta * (1 - cosi * cosi);
    if (k < 0)
        return false;
    else {
        refracted = (eta * I + (eta * cosi - sqrtf(k)) * n).normalized();
        return true;
    }

   /* vec3 uv = (I).normalized();
    float dt = uv.dot(N);
    float d = 1.0 - eta * eta * (1 - dt * dt);
    if (d > 0) {
        refracted = eta * (uv - N * dt) - N * sqrt(d);
        return true;
    }
    else {
        return false;
    }*/


}

//比较函数用以进行bvh划分
bool cmpx(Shape* t1,  Shape* t2) {
    return t1->center.x() < t2->center.x();
}
bool cmpy( Shape* t1,  Shape* t2) {
    return t1->center.y() < t2->center.y();
}
bool cmpz( Shape* t1, Shape* t2) {
    return t1->center.z() < t2->center.z();
}

class BVHNode {
public:
    BVHNode* left = NULL;       // 左右子树索引
    BVHNode* right = NULL;
    int n, index;               // 叶子节点信息               
    vec3 AA, BB;                // 碰撞盒
};
// 构建 BVH
BVHNode* buildBVH(std::vector<Shape*>& shapes, int l, int r, int n) {
    if (l > r) return 0;

    BVHNode* node = new BVHNode();
    node->AA = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
    node->BB = vec3(FLT_MIN, FLT_MIN, FLT_MIN);

    // 计算 AABB
    for (int i = l; i <= r; i++) {
        // 最小点 AA
        vec3 minaa = shapes[i]->minAA();
        node->AA.x() = min(node->AA.x(), minaa.x());
        node->AA.y() = min(node->AA.y(), minaa.y());
        node->AA.z() = min(node->AA.z(), minaa.z());
        // 最大点 BB
        vec3 maxbb = shapes[i]->maxBB();
        node->BB.x() = max(node->BB.x(), maxbb.x());
        node->BB.y() = max(node->BB.y(), maxbb.y());
        node->BB.z() = max(node->BB.z(), maxbb.z());
    }

    // 不多于 n 个三角形 返回叶子节点
    if ((r - l + 1) <= n) {
        node->n = r - l + 1;
        node->index = l;
        return node;
    }

    float Cost = FLT_MAX;
    int Axis = 0;
    int Split = (l + r) / 2;
    for (int axis = 0; axis < 3; axis++) {
        // 分别按 x，y，z 轴排序
        if (axis == 0) std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpx);
        if (axis == 1) std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpy);
        if (axis == 2) std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpz);

        // leftMax[i]: [l, i] 中最大的 xyz 值
        // leftMin[i]: [l, i] 中最小的 xyz 值
        std::vector<vec3> leftMax(r - l + 1, vec3(FLT_MIN, FLT_MIN, FLT_MIN));
        std::vector<vec3> leftMin(r - l + 1, vec3(FLT_MAX, FLT_MAX, FLT_MAX));
        // 计算前缀 注意 i-l 以对齐到下标 0
        for (int i = l; i <= r; i++) {
            int bias = (i == l) ? 0 : 1;  // 第一个元素特殊处理
            vec3 minaa = shapes[i]->minAA();
            vec3 maxbb = shapes[i]->maxBB();

            leftMax[i - l].x() = max(leftMax[i - l - bias].x(), maxbb.x());
            leftMax[i - l].y() = max(leftMax[i - l - bias].y(), maxbb.y());
            leftMax[i - l].z() = max(leftMax[i - l - bias].z(), maxbb.z());

            leftMin[i - l].x() = min(leftMin[i - l - bias].x(), minaa.x());
            leftMin[i - l].y() = min(leftMin[i - l - bias].y(), minaa.y());
            leftMin[i - l].z() = min(leftMin[i - l - bias].z(), minaa.z());
        }

        // rightMax[i]: [i, r] 中最大的 xyz 值
        // rightMin[i]: [i, r] 中最小的 xyz 值
        std::vector<vec3> rightMax(r - l + 1, vec3(FLT_MIN, FLT_MIN, FLT_MIN));
        std::vector<vec3> rightMin(r - l + 1, vec3(FLT_MAX, FLT_MAX, FLT_MAX));
        // 计算后缀 注意 i-l 以对齐到下标 0
        for (int i = r; i >= l; i--) {
            vec3 minaa = shapes[i]->minAA();
            vec3 maxbb = shapes[i]->maxBB();
            int bias = (i == r) ? 0 : 1;  // 第一个元素特殊处理

            rightMax[i - l].x() = max(rightMax[i - l + bias].x(), maxbb.x());
            rightMax[i - l].y() = max(rightMax[i - l + bias].y(), maxbb.y());
            rightMax[i - l].z() = max(rightMax[i - l + bias].z(), maxbb.z());

            rightMin[i - l].x() = min(rightMin[i - l + bias].x(), minaa.x());
            rightMin[i - l].y() = min(rightMin[i - l + bias].y(), minaa.y());
            rightMin[i - l].z() = min(rightMin[i - l + bias].z(), minaa.z());
        }

        // 遍历寻找分割
        float cost = FLT_MAX;
        int split = l;
        for (int i = l; i <= r - 1; i++) {
            float lenx, leny, lenz;
            // 左侧 [l, i]
            vec3 leftAA = leftMin[i - l];
            vec3 leftBB = leftMax[i - l];
            lenx = leftBB.x() - leftAA.x();
            leny = leftBB.y() - leftAA.y();
            lenz = leftBB.z() - leftAA.z();
            float leftS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float leftCost = leftS * (i - l + 1);

            // 右侧 [i+1, r]
            vec3 rightAA = rightMin[i + 1 - l];
            vec3 rightBB = rightMax[i + 1 - l];
            lenx = rightBB.x() - rightAA.x();
            leny = rightBB.y() - rightAA.y();
            lenz = rightBB.z() - rightAA.z();
            float rightS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float rightCost = rightS * (r - i);

            // 记录每个分割的最小答案
            float totalCost = leftCost + rightCost;
            if (totalCost < cost) {
                cost = totalCost;
                split = i;
            }
        }
        // 记录每个轴的最佳答案
        if (cost < Cost) {
            Cost = cost;
            Axis = axis;
            Split = split;
        }
    }

    // 按最佳轴分割
    if (Axis == 0) std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpx);
    if (Axis == 1) std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpy);
    if (Axis == 2) std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpz);

    // 递归
    node->left = buildBVH(shapes, l, Split, n);
    node->right = buildBVH(shapes, Split + 1, r, n);

    return node;

    //// 否则递归建树
    //float lenx = node->BB.x() - node->AA.x();
    //float leny = node->BB.y() - node->AA.y();
    //float lenz = node->BB.z() - node->AA.z();
    //// 按 x 划分
    //if (lenx >= leny && lenx >= lenz)
    //    std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpx);
    //// 按 y 划分
    //if (leny >= lenx && leny >= lenz)
    //    std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpy);
    //// 按 z 划分
    //if (lenz >= lenx && lenz >= leny)
    //    std::sort(shapes.begin() + l, shapes.begin() + r + 1, cmpz);

    //// 递归
    //int mid = (l + r) / 2;
    //node->left = buildBVH(shapes, l, mid, n);
    //node->right = buildBVH(shapes, mid + 1, r, n);

    //return node;
}
float hitAABB(Ray r, vec3 AA, vec3 BB) {//与AABB包围盒的相交
    // 1.0 / direction
   
    vec3 in = (BB - r.startPoint).array() / r.direction.array();
    vec3 out = (AA - r.startPoint).array() / r.direction.array();

    vec3 tmax(max(in.x(),out.x()),max(in.y(),out.y()),max(in.z(),out.z()));// = max(in, out);
    vec3 tmin(min(in.x(), out.x()),min(in.y(), out.y()), min(in.z(), out.z())); //= min(in, out);

    float t1 = min(tmax.x(), min(tmax.y(), tmax.z()));
    float t0 = max(tmin.x(), max(tmin.y(), tmin.z()));

    return (t1 >= t0) ? ((t0 > 0.0) ? (t0) : (t1)) : (-1);
}
//和AABBtree子节点的求交
HitResult hitTriangleArray(Ray ray, std::vector<Shape*>& shapes, int l, int r) {
    HitResult res, tepres;
    res.distance = FLT_MAX;
    for (int i = l; i <= r; i++) {
        tepres =shapes[i]->intersect(ray);
        if (tepres.isHit && tepres.distance < res.distance) {
            res = tepres;
        }
    }
    return res;
}
HitResult hitBVH(Ray ray, std::vector<Shape*>& shapes, BVHNode* root) {
    if (root == NULL) 
        return HitResult();

    // 是叶子 
    if (root->n > 0) {
        return hitTriangleArray(ray, shapes, root->index, root->n + root->index - 1);
    }

    // 和左右子树 AABB 求交
    float d1 = FLT_MAX, d2 = FLT_MAX;
    if (root->left) d1 = hitAABB(ray, root->left->AA, root->left->BB);
    if (root->right) d2 = hitAABB(ray, root->right->AA, root->right->BB);

    // 递归结果
    HitResult r1, r2;
    if (d1 > 0) r1 = hitBVH(ray, shapes, root->left);
    if (d2 > 0) r2 = hitBVH(ray, shapes, root->right);

    return r1.distance < r2.distance ? (r1.isHit?r1:r2) : (r2.isHit?r2:r1);
}