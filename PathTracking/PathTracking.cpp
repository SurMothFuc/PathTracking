// PathTracking.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "common.h"
#include "Secne.h"



void imshow(double* SRC,Scene scene);
HitResult shoot(Scene scene, Ray ray);
vec3 pathTracing(Scene scene, Ray ray,int depth);
void readObj(std::string filepath, Scene& scene);

#define vec3 Vector3f

int main()
{
    // 采样次数
    const int SAMPLE =1;
    // 每次采样的亮度
    const double BRIGHTNESS = (2.0f * 3.1415926f) * (1.0f / double(SAMPLE));//因为是在半球上的均匀采样


    Scene scene;


    readObj("./bunny.obj", scene);


    // 光源
   /* Triangle l1 = Triangle(vec3(0.4, 0.99, 0.4), vec3(-0.4, 0.99, -0.4), vec3(-0.4, 0.99, 0.4), vec3(1.0, 1.0, 1.0));
    Triangle l2 = Triangle(vec3(0.4, 0.99, 0.4), vec3(0.4, 0.99, -0.4), vec3(-0.4, 0.99, -0.4), vec3(1.0, 1.0, 1.0));
    l1.material.isEmissive = true;
    l2.material.isEmissive = true;
    scene.addShape(&l1);
    scene.addShape(&l2);*/
    

    //under
    Triangle* tep;
    tep = new Triangle(vec3(1, -1, 1), vec3(-1, -1, -1), vec3(-1, -1, 1), vec3(1.0, 1.0, 1.0));
    scene.addShape(tep);
    tep = new Triangle(vec3(1, -1, 1), vec3(1, -1, -1), vec3(-1, -1, -1), vec3(1.0, 1.0, 1.0));
    scene.addShape(tep);
    //// top
   
    //tep = new Triangle(vec3(1, 1, 1), vec3(-1, 1, 1), vec3(-1, 1, -1), vec3(1.0, 1.0, 1.0));
    //scene.addShape(tep);
    //tep = new Triangle(vec3(1, 1, 1), vec3(-1, 1, -1), vec3(1, 1, -1), vec3(1.0, 1.0, 1.0));
    //scene.addShape(tep);
    // back
    scene.addShape(new Triangle(vec3(1, -1, -1), vec3(-1, 1, -1), vec3(-1, -1, -1), vec3(1.0, 1.0, 1.0)));
    scene.addShape(new Triangle(vec3(1, -1, -1), vec3(1, 1, -1), vec3(-1, 1, -1), vec3(1.0, 1.0, 1.0)));
    //// left
    //scene.addShape(new Triangle(vec3(-1, -1, -1), vec3(-1, 1, 1), vec3(-1, -1, 1), vec3(1.0, 0.0, 0.0)));
    //scene.addShape(new Triangle(vec3(-1, -1, -1), vec3(-1, 1, -1), vec3(-1, 1, 1), vec3(1.0, 0.0, 0.0)));
    //// right
    //scene.addShape(new Triangle(vec3(1, 1, 1), vec3(1, -1, -1), vec3(1, -1, 1), vec3(0.0, 1.0, 0.0)));
    //scene.addShape(new Triangle(vec3(1, -1, -1), vec3(1, 1, 1), vec3(1, 1, -1), vec3(0.0, 1.0, 0.0)));


    double* image = (double*)malloc(sizeof(double) * scene.WIDTH * scene.HEIGHT * 3);
    memset(image, 0.0, sizeof(double) * scene.WIDTH * scene.HEIGHT * 3);
    /*   for (int j = 0; j < 720; j++) {
        for (int i = 0; i < 1280; i++) {
            s[j*1280*3+3*i] = 1.0*i/1280;
            s[j * 1280 * 3 + i* 3+1] = 1.0*j/720;
        } 
    }*/
    vec3 EYE(0,0,4.0);



    static omp_lock_t lock;
    omp_init_lock(&lock);


    ProcessBar processbar(SAMPLE, "path tracing:");
    processbar.start();
    int count = 0;
    omp_set_num_threads(8); // 线程个数
    #pragma omp parallel for
    for (int k = 0; k < SAMPLE; k++)
    {
       // cout << k << endl;
        double* p = image;
        for (int i = 0; i < scene.HEIGHT; i++) {
            for (int j = 0; j < scene.WIDTH; j++) {
                double x = 2.0 * double(j) / double(scene.WIDTH) - 1;//映射到-1 1 投影平面的范围内
                double y = 2.0 * double(scene.HEIGHT - i) / double(scene.HEIGHT) - 1;//因为图片存储是从左上角开始 所以这里翻转一下y轴

                vec3 coord(x, y, scene.SCREEN_Z);//在投影平面上的坐标
                vec3 direction = (coord - EYE).normalized();


                Ray ray;
                ray.startPoint = coord;
                ray.direction = direction;
                HitResult res = shoot(scene, ray);

                vec3 color = vec3(0.0, 0.0, 0.0);

                if (res.isHit)
                {
                    // 命中光源直接返回光源颜色
                    if (res.material.isEmissive)
                    {
                        color = res.material.color;
                    }
                    // 命中实体则选择一个随机方向重新发射光线并且进行路径追踪
                    else
                    {
                        // 根据交点处法向量生成交点处反射的随机半球向量
                        Ray randomRay;
                        randomRay.startPoint = res.hitPoint;
                        randomRay.direction = randomDirection(res.material.normal);


                        double r = randf();
                        if (r < res.material.specularRate)  // 镜面反射
                        {
                            randomRay.direction = reflect(ray.direction, res.material.normal).normalized();
                            color = pathTracing(scene, randomRay, 0) * BRIGHTNESS;
                        }
                        // 颜色积累
                        else {
                            vec3 srcColor = res.material.color;
                            vec3 ptColor = pathTracing(scene, randomRay, 0);
                            color = ptColor.array() * srcColor.array() * BRIGHTNESS;    // 和原颜色混合
                        }

                    }
                }

               
                omp_set_lock(&lock);
                *p += color.x(); p++;  // R 通道
                *p += color.y(); p++;  // G 通道
                *p += color.z(); p++;  // B 通道
                omp_unset_lock(&lock); //释放互斥器

            }
        }
        omp_set_lock(&lock);
        count++;
        processbar.update(count + 1);
        omp_unset_lock(&lock); //释放互斥器
    }

    processbar.end();
    imshow(image,scene);
    omp_destroy_lock(&lock);
}


void readObj(std::string filepath,Scene& scene) {
    // 打开文件流
    std::vector<vec3> vertices;
    std::ifstream fin(filepath);
    std::string line;
    if (!fin.is_open()) {
        std::cout << "文件 " << filepath << " 打开失败" << std::endl;
        exit(-1);
    }

    // 增量读取
    int offset = vertices.size();
    int count = 0;
    // 按行读取
    while (std::getline(fin, line)) {
        std::istringstream sin(line);   // 以一行的数据作为 string stream 解析并且读取
        std::string type;
        float x, y, z;
        int v0, v1, v2;

        // 读取obj文件
        sin >> type;
        if (type == "v") {
            sin >> x >> y >> z;
            vec3 point = 5 * vec3(x, y, z);
            point.z() += 1.0;
            point.y() -= 0.9;
            point.x() += 0.09;
            vertices.push_back(point);
        }
        if (type == "f") {
            sin >> v0 >> v1 >> v2;
            Triangle* tep;
            tep = new Triangle(vertices[v0-1], vertices[v1 - 1], vertices[v2 - 1], vec3(1.0, 1.0, 1.0));

            tep->material.isEmissive = true;
            scene.addShape(tep);

        }
    }
}





vec3 pathTracing(Scene scene, Ray ray, int depth)
{
    if (depth > 0) 
        return vec3(0.0,0.0,0.0);

    HitResult res = shoot(scene, ray);

    if (!res.isHit) 
        return vec3(0.0, 0.0, 0.0); // 未命中

    // 如果发光则返回颜色
    if (res.material.isEmissive) 
        return res.material.color;


    // 俄罗斯转盘
    double r = randf();
    float P = 0.8;
    if (r > P) 
        return vec3(0.0, 0.0, 0.0);

    Ray randomRay;
    randomRay.startPoint = res.hitPoint;
    randomRay.direction = randomDirection(res.material.normal);

    float cosine = (-ray.direction).dot(res.material.normal);
    vec3 color(0, 0, 0);

    r = randf();
    if (r < res.material.specularRate)  // 镜面反射
    {
        randomRay.direction = reflect(ray.direction, res.material.normal).normalized();
        color = pathTracing(scene, randomRay, depth + 1) * cosine;
    }
    else    // 漫反射
    {
        vec3 srcColor = res.material.color;
        vec3 ptColor = pathTracing(scene, randomRay, depth + 1) * cosine;
        color = ptColor.array() * srcColor.array();    // 和原颜色混合
    }
    return color / P;
}

HitResult shoot(Scene scene, Ray ray)
{
    HitResult res, r;
    res.distance = FLT_MAX; // inf

    // 遍历所有图形，求最近交点
    for (auto& shape : scene.shapes)
    {
        r = shape->intersect(ray);
        if (r.isHit && r.distance < res.distance)
            res = r;  // 记录距离最近的求交结果
    }

    return res;
}
void imshow(double* SRC,Scene scene)
{

    unsigned char* image = new unsigned char[scene.WIDTH * scene.HEIGHT * 3];// 图像buffer
    unsigned char* p = image;
    double* S = SRC;    // 源数据

    FILE* fp;
    fopen_s(&fp, "image.png", "wb");

    for (int i = 0; i < scene.HEIGHT; i++)
    {
        for (int j = 0; j < scene.WIDTH; j++)
        {
            //进行下伽马校正 不然暗的地方会太暗
            *p++ = (unsigned char)Clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0);  // R 通道
            *p++ = (unsigned char)Clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0);  // G 通道
            *p++ = (unsigned char)Clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0);  // B 通道
        }
    }

    svpng(fp, scene.WIDTH, scene.HEIGHT, image, 0);
}