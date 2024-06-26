#include <iostream>
#include <vector>
#include <cmath> // For math functions
#include <ctime> // For time
#include <cstdlib> // For srand and rand
#include <unordered_map>
#include <cmath>
#include <unordered_set> 
#include <array>
#include <random>
#include <map>
#include <set>
#include "libxl.h" // 包含 LibXL 头文件
using namespace libxl;
using namespace std;

#define pi acos(-1)

// 计算距离
double dis(std::vector<double>& a, std::vector<double>& b) {
    double sum2 = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum2 += pow((a[i] - b[i]), 2);
    }
    return sqrt(sum2);
}

// 读取 Excel 文件并返回二维字符串向量
std::vector<std::vector<double>> read_data(const wchar_t* path) {
    Book* book = xlCreateBook(); // 创建一个 Excel 文档对象
    std::vector<std::vector<double>> data;

    if (book->load(path)) { // 加载 Excel 文件
        Sheet* sheet = book->getSheet(0); // 获取第一个工作表
        if (sheet) {
            //int rowCount = sheet->lastRow(); // 获取行数
            //int colCount = sheet->lastCol(); // 获取列数

            int rowCount = sheet->lastRow(); // 获取行数
            int colCount = sheet->lastCol(); // 获取列数
            for (int i = 0; i <= rowCount; ++i) {
                std::vector<double> rowData;
                for (int j = 0; j <= colCount; ++j) {
                    double value = sheet->readNum(i, j); // 读取单元格数据
                    std::cout << "输出第"<<i<<"行第"<<j<<"列" << value << "\t"; // 输出数据
                    rowData.push_back(value);
                }
                data.push_back(rowData);
            }
        }
        else {
            std::cerr << "Error accessing sheet!" << std::endl;
        }
    }
    else {
        std::cerr << "Error loading Excel file!" << std::endl;
    }

    book->release(); // 释放 Excel 文档对象
    return data;
}

//打印数组
/* 使用方式
    int rows = sizeof(dan_list) / sizeof(dan_list[0]);
    int cols = sizeof(dan_list[0]) / sizeof(dan_list[0][0]);

    printArray(dan_list, rows, cols);  */
void printArray(double arr[][4], int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}



//定义Target类

class Target {
public:
    Target(int number) : number(number), location({ 80., 20., 0. }), theta(pi), velocity(0.016), value(0),
        health(1), bullet_num(0), bullet_real(0), theta_head(pi / 6), theta_side(pi / 6),
        theta_end(pi / 9), theta_domain0({ -pi + theta_end, -pi / 2 - theta_side }),
        theta_domain1({ -pi / 2 + theta_side, -theta_head }), theta_domain2({ theta_head, pi / 2 - theta_side }),
        theta_domain3({ pi / 2 + theta_side, pi - theta_end }) {}

    void update_target(std::vector<double> location, int bullet_num, int health) {
        this->location = location;
        this->health = health;
        this->bullet_num = bullet_num;
    }

    void moving(double dt) {
        this->location[0] -= this->velocity * dt;
    }

public:
    int number;
    std::vector<double> location;
    double theta;
    double velocity;
    int value;
    int health;
    int bullet_num;
    int bullet_real;
    double theta_head;
    double theta_side;
    double theta_end;
    std::vector<double> theta_domain0;
    std::vector<double> theta_domain1;
    std::vector<double> theta_domain2;
    std::vector<double> theta_domain3;
};

// 定义dubins距离 和经纬度转换
double rad_normol(double theta) {
    if (theta > pi) {
        theta -= 2 * pi;
    }
    else if (theta < -pi) {
        theta += 2 * pi;
    }
    return theta;
}

double dis_lon_lat(double lon1, double lat1, double lon2, double lat2) {
    double R = 6371;  // 地球半径（千米）
    // 度转弧度
    double lon1rad = lon1 * pi / 180;
    double lat1rad = lat1 * pi / 180;
    double lon2rad = lon2 * pi / 180;
    double lat2rad = lat2 * pi / 180;
    double d = R * acos(cos(lat1rad) * cos(lat2rad) * cos(lon1rad - lon2rad) + sin(lat1rad) * sin(lat2rad));
    // 百公里误差数米，千公里误差数十米
    return d;
}

double mod(double a, double b) {
    return fmod(a, b);
}

std::array<double, 4> lsl(double alpha, double beta, double d) {
    double p_squared = 2 + (d * d) - (2 * cos(alpha - beta)) + (2 * d * (sin(alpha) - sin(beta)));
    if (p_squared < 0) {
        return { INFINITY, INFINITY, INFINITY, INFINITY };
    }

    double tmp0 = d + sin(alpha) - sin(beta);
    double tmp1 = atan2((cos(beta) - cos(alpha)), tmp0);
    double t = mod((-alpha + tmp1), 2 * pi);
    double p = sqrt(p_squared);
    double q = mod((beta - tmp1), 2 * pi);
    double length = t + p + q;

    return { length, t, p, q };
}

std::array<double, 4> lsr(double alpha, double beta, double d) {
    double p_squared = -2 + d * d + 2 * cos(alpha - beta) + 2 * d * (sin(alpha) + sin(beta));
    if (p_squared < 0) {
        return { INFINITY, INFINITY, INFINITY, INFINITY };
    }
    double p = sqrt(p_squared);
    double tmp2 = atan2((-cos(alpha) - cos(beta)), (d + sin(alpha) + sin(beta))) - atan2(-2.0, p);
    double t = mod((-alpha + tmp2), 2 * pi);
    double q = mod((-mod(beta, 2 * pi) + tmp2), 2 * pi);
    double length = t + p + q;

    return { length, t, p, q };
}

std::array<double, 4> rsl(double alpha, double beta, double d) {
    double p_squared = (d * d) - 2 + (2 * cos(alpha - beta)) - (2 * d * (sin(alpha) + sin(beta)));
    if (p_squared < 0) {
        return { INFINITY, INFINITY, INFINITY, INFINITY };
    }
    double p = sqrt(p_squared);
    double tmp2 = atan2((cos(alpha) + cos(beta)), (d - sin(alpha) - sin(beta))) - atan2(2.0, p);
    double t = mod((alpha - tmp2), 2 * pi);
    double q = mod((beta - tmp2), 2 * pi);
    double length = t + p + q;

    return { length, t, p, q };
}

std::array<double, 4> rsr(double alpha, double beta, double d) {
    double p_squared = 2 + (d * d) - (2 * cos(alpha - beta)) + (2 * d * (sin(beta) - sin(alpha)));
    if (p_squared < 0) {
        return { INFINITY, INFINITY, INFINITY, INFINITY };
    }
    double tmp0 = d - sin(alpha) + sin(beta);
    double tmp1 = atan2((cos(alpha) - cos(beta)), tmp0);
    double t = mod((alpha - tmp1), 2 * pi);
    double p = sqrt(p_squared);
    double q = mod((-beta + tmp1), 2 * pi);
    double length = t + p + q;

    return { length, t, p, q };
}

std::array<double, 4> rlr(double alpha, double beta, double d) {
    double tmp_rlr = (6. - d * d + 2 * cos(alpha - beta) + 2 * d * (sin(alpha) - sin(beta))) / 8.;
    if (std::abs(tmp_rlr) > 1) {
        return { INFINITY, INFINITY, INFINITY, INFINITY };
    }
    double p = mod((2 * pi - acos(tmp_rlr)), 2 * pi);
    double t = mod((alpha - atan2(cos(alpha) - cos(beta), d - sin(alpha) + sin(beta)) + mod(p / 2, 2 * pi)), 2 * pi);
    double q = mod((alpha - beta - t + mod(p, 2 * pi)), 2 * pi);
    double length = t + p + q;

    return { length, t, p, q };
}

std::array<double, 4> lrl(double alpha, double beta, double d) {
    double tmp_lrl = (6. - d * d + 2 * cos(alpha - beta) + 2 * d * (-sin(alpha) + sin(beta))) / 8.;
    if (std::abs(tmp_lrl) > 1) {
        return { INFINITY, INFINITY, INFINITY, INFINITY };
    }
    double p = mod((2 * pi - acos(tmp_lrl)), 2 * pi);
    double t = mod((-alpha - atan2(cos(alpha) - cos(beta), d + sin(alpha) - sin(beta)) + p / 2), 2 * pi);
    double q = mod((mod(beta, 2 * pi) - alpha - t + mod(p, 2 * pi)), 2 * pi);
    double length = t + p + q;

    return { length, t, p, q };
}

double dubins_min_len( std::array<double, 3>& p1, double theta1,  std::vector<double>& p2, double theta2, double r) {
    // 经纬度距离转化，建立相对平面直角坐标系
    // double dx = dis_lon_lat(p1[0], p2[1], p2[0], p2[1]);
    // double dy = dis_lon_lat(p2[0], p1[1], p2[0], p2[1]);

    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double d = sqrt(dx * dx + dy * dy) / r;

    double theta = mod(atan2(dy, dx), 2 * pi);
    double alpha = mod((theta1 - theta), 2 * pi);
    double beta = mod((theta2 - theta), 2 * pi);

    std::array<double, 4> length = lsl(alpha, beta, d);
    std::array<double, 4> tmp = lsr(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
    }
    tmp = rsl(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
    }

    tmp = rsr(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
    }

    tmp = rlr(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
    }

    tmp = lrl(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
    }

    return r * length[0];
}

vector<double> dubins_segment(double seg_param, vector<double> seg_init, char seg_type) {
    vector<double> seg_end(3, 0);

    if (seg_type == 'L') {
        seg_end[0] = seg_init[0] + sin(seg_init[2] + seg_param) - sin(seg_init[2]);
        seg_end[1] = seg_init[1] - cos(seg_init[2] + seg_param) + cos(seg_init[2]);
        seg_end[2] = seg_init[2] + seg_param;
    }
    else if (seg_type == 'R') {
        seg_end[0] = seg_init[0] - sin(seg_init[2] - seg_param) + sin(seg_init[2]);
        seg_end[1] = seg_init[1] + cos(seg_init[2] - seg_param) - cos(seg_init[2]);
        seg_end[2] = seg_init[2] - seg_param;
    }
    else if (seg_type == 'S') {
        seg_end[0] = seg_init[0] + cos(seg_init[2]) * seg_param;
        seg_end[1] = seg_init[1] + sin(seg_init[2]) * seg_param;
        seg_end[2] = seg_init[2];
    }

    return seg_end;
}
vector<double> dubins_end_point(std::array<double, 3>&  p1, double theta1, vector<double> p2, double theta2, double r, double s) {

    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double d = sqrt(dx * dx + dy * dy) / r;

    double theta = fmod(atan2(dy, dx), 2 * pi);
    double alpha = fmod((theta1 - theta), 2 * pi);
    double beta = fmod((theta2 - theta), 2 * pi);

    std::array<double, 4> length = lsl(alpha, beta, d);
    char type_d = 'LSL';
    std::array<double, 4> tmp = lsr(alpha, beta, d);

    if (tmp[0] < length[0]) {
        length = tmp;
        type_d = 'LSR';
    }
    tmp = rsl(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
        type_d = 'RSL';
    }
    tmp = rsr(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
        type_d = 'RSR';
    }
    tmp = rlr(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
        type_d = 'RLR';
    }
    tmp = lrl(alpha, beta, d);
    if (tmp[0] < length[0]) {
        length = tmp;
        type_d = 'LRL';
    }

    double seg_param = s / r;
    vector<double> seg_init = { 0, 0, theta1 };
    if (seg_param <= length[1]) {
        vector<double> seg_end = dubins_segment(seg_param, seg_init, type_d);
        double x = p1[0] + seg_end[0] * r;
        double y = p1[1] + seg_end[1] * r;
        theta = fmod(seg_end[2], 2 * pi);
        return { x, y, theta };
    }

    vector<double> mid1 = dubins_segment(length[1], seg_init, type_d);
    if (seg_param <= length[1] + length[2]) {
        vector<double> seg_end = dubins_segment(seg_param - length[1], mid1, type_d);
        double x = p1[0] + seg_end[0] * r;
        double y = p1[1] + seg_end[1] * r;
        theta = fmod(seg_end[2], 2 * pi);
        return { x, y, theta };
    }

    vector<double> mid2 = dubins_segment(length[2], mid1, type_d);
    vector<double> seg_end = dubins_segment(seg_param - length[1] - length[2], mid2, type_d);
    double x = p1[0] + seg_end[0] * r;
    double y = p1[1] + seg_end[1] * r;
    theta = fmod(seg_end[2], 2 * pi);
    return { x, y, theta };
}

//定义消息类
class Message {
public:
    int receiver;  // 接收方编号
    int sender;    // 发送方编号
    int send_time; // 当前时间
    int receive_time;  // 在环境中传播后，接收到的时间
    double delay;  // 环境时延
    std::vector<int> text;  // 信息内容

public:
    Message(int _receiver, int _sender, int time, std::vector<int>_text) :
        receiver(_receiver), sender(_sender), send_time(time), text(_text) {
        // 计算环境时延
        // double delay = 3 + 3 * sin(0.05 * time + 1) + 2 * sin(1.22 * time + .5 * time + 0.2) + 1.5 * sin(2.2 * time + 1) + cos(1.8 * time + 1) ** 2;
        // if (delay < 1)
        //     delay = 1;
        delay = 1.0;
        receive_time = time + static_cast<int>(delay);
    }

    // 其他成员函数...
};
//定义导弹类
class Missile {
public:
    Missile(int _number) : number(_number) {}
public:
    int number;
    std::array<double, 3> location;
    double s_range = 400.0;  // 剩余航程 km
    double velocity = 0.200;  // 速度 km/s
    double radius = 20;  // 转弯半径 km
    double theta = 0;  // 航向角 rad
    double hit_theta = 0;  // 打击入射角 rad
    int value = 100;  // 价值
    double commu_distance = 20;  // 通信距离 km
    double commu_prob = 1;  // 通信质量 0~1之间
    int health = 1;  // 当前状态 0 损毁， 1 健康
    int hit = 0;  // 是否击中目标， 0 未击中， 1 击中
    double angle_min = pi / 60;  // 相邻弹最小打击角度间隔，3度
    Target task1 = Target(0);  // 执行任务信息
    Target task2 = Target(0);  // 携带任务信息
    double advantage_val = 0;  // 打击优势值
    double hit_distance = 0;  // 距离目标剩余打击距离
   // std::map<double, std::map<double, double>> hit_theta_dict; // 用于存放打击同一任务成员的打击角度值,key 1 为任务编号,key 2为导弹编号,值为打击角度
    std::map<int, std::map<int, std::vector<double>>> hit_theta_dict;
    std::vector<std::pair<int, std::vector<double>>> hit_theta_map;
    //std::vector<std::pair<double, double>> hit_theta_map; // 用于存放打击同一任务成员的打击角度值


    //std::unordered_map<int, double> hit_theta_dict; 
    //std::vector<double> hit_theta_map;  // 用于存放打击同一任务成员的打击角度值
    std::array<double, 3> hit_theta_left = { -10, 0, 0 };  // 用于存放打击同一任务成员的打击角度值紧临的小于值,存储格式[角度，导弹编号，任务编号]
    std::array<double, 3> hit_theta_right = { 10, 0, 0 };  // 用于存放打击同一任务成员的打击角度值紧临的大于值,存储格式[角度，导弹编号，任务编号]
    std::unordered_map<int, double> task_value_dict;  // 用于存放打击任务成员的剩余航程优势函数值，key 1 为任务编号,key 2为导弹编号
    double task_value = 0;  // 执行任务的优势值，剩余航程
    int phase = 0;  // 通信阶段 0 为第一次通信，1为发送决策内容，2为决策确认
    std::vector<int> message;  // 报文内容
    std::unordered_map<int, int> connect_buffer;  // 连接缓存器，存放保持连接的飞行器，用飞行器编号作为索引
    std::unordered_map<int, double> state_time;  // 最新接受到飞行器消息的时间，用飞行器编号作为索引，用于判断飞行器是否存活
    std::unordered_map<int, double> connect_time;  // 连接时间，飞行器状态最新更新时间，用飞行器编号作为索引
    std::vector<Target> task_buffer;  // 从盟友回收的任务，或侦察获得的额外任务
    std::unordered_map<int, double> interrupt_buffer;  // 中断缓存器
    std::unordered_map<int, double> hit_buffer;  // 打击角缓存器，用于存在打击相同目标弹群的打击角，键为导弹编号，值为打击角度；按照打击角度进行排序，便于解除冲突
    std::unordered_map<int, int> ally_buffer;  // 协作联盟，用于存放盟友的编号，盟友的信息数据在connect_buffer中，键为盟友编号，值为盟友状态，0为正常，
    // 1 2 3为发送脱离信息次数

    std::vector<double> message_buffer;  // 用于暂放消息
    int buy_stage = 0;  // 1:广播任务，2：投标，3：中标结果，4：确认执行
    int buy_flag = 0;  // 标识初始化
    double buy_time = 0;  // 计时器初始化
    int buy_number = 0;  // 存放交易对象编号
    int buy_need = 0;  // 买入需求:1为有需求，0为无需求

    double buy_time_begin = 2;
    double buy_deadline = 10;
    Target buy_target = Target(0);
    std::unordered_map<int, double> buy_bidder_buffer;  // 存放标书

    std::unordered_map<int, Target> buy_task_buffer;  // 存放接收到的任务，安照发送方编号存储
    std::unordered_map<int, std::tuple<double, double, double>> buy_task_value; // 存放接收到的任务优势值，安照发送方编号存储
    std::unordered_set<int> buy_reply_buffer;  // 等待回复
    double buy_reply_timer = 0;  // 等待回复计时器
    int buy_reply_flag = 0;  // 等待计时标志
    int buy_reply_message = 0;  // 需要重发的信息

    int sell_stage = 0;
    double sell_time = 0;
    int sell_flag = 0;  // 0:买，1：卖
    int sell_need = 0;
    int sell_num = 0;  // 交易对象编号
    int sell_bullet_num = 0;
    std::unordered_set<int> sell_reply_buffer;  // 等待回复
    double sell_reply_timer = 0;  // 等待回复计时器
    int sell_reply_flag = 0;  // 等待计时标志
    int sell_reply_message = 0;  // 需要重发的信息

    int exchange_number = 0;
    int exchange_stage = 0;  // 1：请求交换，2：回复请求，3：确认执行
    double exchange_time = 0;
    int exchange_flag = 0;
    double exchange_deadline = 10;
    Target exchange_target = Target(0);
    std::unordered_set<int> exchange_reply_buffer;  // 等待回复
    double exchange_reply_timer = 0;  // 等待回复计时器
    int exchange_reply_flag = 0;  // 等待计时标志
    int exchange_reply_message = 0;  // 需要重发的信息
    std::vector<Target> target_buffer;  // 暂存用于交换的目标
    std::vector<std::string> exchange_text;  // 用于存放交换的信息

    // 重置任务相关数据
    void task_reset() {
        hit_theta_left = { -10, 0, 0 }; // 用于存放打击同一任务成员的打击角度值紧临的小于值,存储格式[角度，导弹编号，任务编号]
        hit_theta_right = { 10, 0, 0 }; // 用于存放打击同一任务成员的打击角度值紧临的大于值,存储格式[角度，导弹编号，任务编号]
    }

    // 重置购买相关数据
    void buy_reset() {
        buy_stage = 0; // 1:广播任务，2：投标，3：中标结果，4：确认执行
        buy_flag = 0; // 标识初始化
        buy_time = 0; // 计时器初始化
        buy_number = 0; // 存放交易对象编号
        buy_need = 0; // 买入需求:1为有需求，0为无需求
        buy_target.number = 0; // 交易任务初始化
    }

    // 重置出售相关数据
    void sell_reset() {
        sell_stage = 0;
        sell_time = 0;
        sell_flag = 0; // 0:买，1：卖
    }

    void update_missile( std::array<double, 3>& _location, double _s_range, double _commu_prob, int _health,  Target& _task1,  Target& _task2) {
        location = _location;
        s_range = _s_range;
        commu_prob = _commu_prob;
        health = _health;
        task1 = _task1;
        task2 = _task2;
    }

    double advantage(Target& ship, double hit_theta = 0.0, double distance = 0.0) {
        // 距离优势
        if (hit_theta == 0.0 && distance == 0.0) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);

            hit_theta = ship.theta_domain0[0] + (ship.theta_domain0[1] - ship.theta_domain0[0]) * dis(gen);
            double hit_theta_abs = rad_normol(ship.theta + hit_theta);  // 绝对角度
            distance = dubins_min_len(location, theta, ship.location, hit_theta_abs, radius);

            double hit_theta_tmp = ship.theta_domain1[0] + (ship.theta_domain1[1] - ship.theta_domain1[0]) * dis(gen);
            double hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp);
            double distance_tmp = dubins_min_len(location, theta, ship.location, hit_theta_abs_tmp, radius);
            if (distance_tmp < distance) {
                hit_theta = hit_theta_tmp;
                distance = distance_tmp;
            }
            hit_theta_tmp = ship.theta_domain2[0] + (ship.theta_domain2[1] - ship.theta_domain2[0]) * dis(gen);
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp);
            distance_tmp = dubins_min_len(location, theta, ship.location, hit_theta_abs_tmp, radius);
            if (distance_tmp < distance) {
                hit_theta = hit_theta_tmp;
                distance = distance_tmp;
            }
            hit_theta_tmp = ship.theta_domain3[0] + (ship.theta_domain3[1] - ship.theta_domain3[0]) * dis(gen);
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp);
            distance_tmp = dubins_min_len(location, theta, ship.location, hit_theta_abs_tmp, radius);
            if (distance_tmp < distance) {
                hit_theta = hit_theta_tmp;
                distance = distance_tmp;
            }
        }
        else if (distance == 0.0) {
            double hit_theta_abs = rad_normol(ship.theta + hit_theta);  // 绝对角度
            distance = dubins_min_len(location, theta, ship.location, hit_theta_abs, radius);
        }

        if (distance > s_range) {
            return 0;
        }
        else if (s_range < 50) {
            return 1 + 500 / distance * (ship.value / value);
        }
        else {
            return 1 + distance / s_range + s_range / distance * (ship.value / value);
        }
    }



    std::tuple<double, double, double> decision_making(Target& ship, double hit_theta = -1) {
        // 距离优势
        double distance;
        double vd;
        if (hit_theta < 0) {
            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            double hit_theta = ship.theta_domain0[0] + (ship.theta_domain0[1] - ship.theta_domain0[0]) * distribution(generator);
            double hit_theta_abs = rad_normol(ship.theta + hit_theta);  // 绝对角度
             distance = dubins_min_len(location, theta, ship.location, hit_theta_abs, radius);

            double hit_theta_tmp = ship.theta_domain1[0] + (ship.theta_domain1[1] - ship.theta_domain1[0]) * distribution(generator);
            double hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp);
            double distance_tmp = dubins_min_len(location, theta, ship.location, hit_theta_abs_tmp, radius);
            if (distance_tmp < distance) {
                hit_theta = hit_theta_tmp;
                distance = distance_tmp;
            }

            hit_theta_tmp = ship.theta_domain2[0] + (ship.theta_domain2[1] - ship.theta_domain2[0]) * distribution(generator);
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp);
            distance_tmp = dubins_min_len(location, theta, ship.location, hit_theta_abs_tmp, radius);
            if (distance_tmp < distance) {
                hit_theta = hit_theta_tmp;
                distance = distance_tmp;
            }

            hit_theta_tmp = ship.theta_domain3[0] + (ship.theta_domain3[1] - ship.theta_domain3[0]) * distribution(generator);
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp);
            distance_tmp = dubins_min_len(location, theta, ship.location, hit_theta_abs_tmp, radius);
            if (distance_tmp < distance) {
                hit_theta = hit_theta_tmp;
                distance = distance_tmp;
            }
        }
        else {
            double hit_theta_abs = rad_normol(ship.theta + hit_theta);  // 绝对角度
            distance = dubins_min_len(location, theta, ship.location, hit_theta_abs, radius);
        }

        if (distance > s_range) {
            return std::make_tuple(0.0, 0.0, 0.0);
        }
        else if (s_range < 50) {
           vd = 1.0 + 500.0 / distance;
        }
        else {
            vd = 1.0 + distance / s_range + s_range / distance;
        }

        // 代价优势
        double vv = ship.value / value;

        return std::make_tuple(vd * vv, hit_theta, distance);
    }

    void moving1(double dt) {
        double s = velocity * dt;

        vector<double> x_y_theta= dubins_end_point(location, theta, task1.location, task1.theta + hit_theta, radius, s);
        location[0] = x_y_theta[0];
        location[1] = x_y_theta[1];
        this->theta = x_y_theta[2];
    }

    void moving2(double dt) {
        if (task1.number > 0) { // 有打击目标
            // 调用 decision_making 方法获取 val、hit_theta 和 distance
            double val, hit_theta, distance;
            std::tie(val, hit_theta, distance) = decision_making(task1);

            // 调用 update_hit_theta_map 方法
            update_hit_theta_map();

            if (distance < hit_distance && check_new_hit_theta(hit_theta)) {
                hit_theta = hit_theta;
                hit_distance = distance;
                // task_reset();
            }

            // 调用 angle_check 方法
            // angle_check();

            // 更新 hit_distance 和 advantage_val
            hit_distance = distance;
            advantage_val = advantage(task1, distance = hit_distance);
        }
    }

    void buyer_age() {  // 买方时效性评估
        buy_time += 1;
        if (buy_time > buy_deadline) {
            buy_flag = 0;
            buy_stage = 0;
            buy_time = 0;
            buy_number = 0;
        }
    }

    void seller_age() {
        sell_time += 1;
        if (sell_time > buy_deadline) {
            sell_flag = 0;
            sell_stage = 0;
            sell_time = 0;
        }
    }

    Message broadcast(int dan2_number, int time) {
        if (sell_stage == 0 && task2.number > 0) {
            sell_bullet_num = task2.bullet_num;
            sell_num = 0;
            sell_flag = true;

            std::vector<int> text = { 1, 1,task2.number }; // Assuming the first two elements are some identifiers

            return Message(dan2_number, number, time, text);
        }
        return Message(0, 0, 0, {});
    }

    bool check_new_hit_theta(double new_hit_theta) {
        std::array<double, 3> hit_theta_left = { -10.0, 0.0, 0.0 };
        std::array<double, 3> hit_theta_right = { 10.0, 0.0, 0.0 };

        for (size_t i = 0; i < hit_theta_map.size(); ++i) {
            if (hit_theta_map[i].second[0] > new_hit_theta) {
                hit_theta_right = { hit_theta_map[i].second[0], static_cast<double>(hit_theta_map[i].first), static_cast<double>(task1.number) };
                if (i > 0) {
                    hit_theta_left = { hit_theta_map[i - 1].second[0], static_cast<double>(hit_theta_map[i - 1].first),  static_cast<double>(task1.number) };
                }
                break;
            }
            if (i == hit_theta_map.size() - 1) {
                hit_theta_left = { hit_theta_map[i].second[0], static_cast<double>(hit_theta_map[i].first),  static_cast<double>(task1.number)};
            }
        }

        if (hit_theta_right[0] - hit_theta_left[0] >= 2 * angle_min) {
            this->hit_theta_left = hit_theta_left;
            this->hit_theta_right = hit_theta_right;
            return true;
        }
        else {
            return false;
        }
    }
    void update_hit_theta_map() {
        auto it = hit_theta_dict.find(task1.number);
        if (it != hit_theta_dict.end()) {
            const auto& dict = it->second;
            hit_theta_map.clear();
            for (const auto& item : dict) {
                hit_theta_map.push_back({ item.first, item.second });
            }
            std::sort(hit_theta_map.begin(), hit_theta_map.end(), [](const auto& a, const auto& b) { return a.second[0] < b.second[0]; });
            for (size_t i = 0; i < hit_theta_map.size(); ++i) {
                if (hit_theta_map[i].second[0] > hit_theta) {
                    hit_theta_right = { hit_theta_map[i].second[0], static_cast<double>(hit_theta_map[i].first), static_cast<double>(task1.number) };
                    if (i > 0) {
                        hit_theta_left = { hit_theta_map[i - 1].second[0], static_cast<double>(hit_theta_map[i - 1].first), static_cast<double>(task1.number) };
                    }
                    else {
                        hit_theta_left = { -10.0, 0, 0 };
                    }
                    return;
                }
            }
            if (hit_theta_map.empty()) {
                hit_theta_left = { 10.0, 0, 0 };
                hit_theta_right = { 10.0, 0, 0 };
            }
            else {
                hit_theta_left = { hit_theta_map.back().second[0], static_cast<double>(hit_theta_map.back().first), static_cast<double>(task1_number) };
                hit_theta_right = { 10.0, 0, 0 };
            }
        }
    }

    Message buy_doing(int time) {
        if (buy_flag == 0) {
            buy_time++;
            if (buy_time >= buy_time_begin) {
                buy_time = 0;
                if (!buy_task_value.empty()) {
                    auto it = std::min_element(buy_task_value.begin(), buy_task_value.end(),
                        [](const auto& a, const auto& b) { return a.second < b.second; });
                    buy_number = it->first;
                    buy_target = buy_task_buffer[buy_number];

                    int task1_flag = 0, task2_flag = 0;
                    double val, hit_theta, hit_distance;
                    std::tie(val, hit_theta, hit_distance) = buy_task_value[buy_number];

                    if (buy_target.number == 0 && val > 0) {
                        task1_flag = buy_target.number;
                        advantage_val = val;
                        this->hit_theta = hit_theta;
                        this->hit_distance = hit_distance;
                        buy_target.bullet_num -= 1;
                    }

                    if (task1_flag || task2_flag) {
                        buy_reply_buffer.insert(buy_number);
                        buy_stage = 2;
                        buy_flag = 1;
                        std::vector<int> text = { 1, 2, task1_flag, task2_flag };
                        return Message(buy_number, number, time, text);
                    }
                    else {
                        return Message(0, 0, 0, std::vector<int>());
                    }
                }
            }
        }
        return Message(0, 0, 0, std::vector<int>());
    }

};



int main() {
    // 定义 con_para 字典
    std::unordered_map<std::string, double> con_para = {
        {"com_pro", 1},
        {"shoot_pro", 1}
    };

    // 初始化数据
    int time_slot = 20;
    int time_end = 3001;
    int dt = 1;  // 飞行器采样频率


    //（1）读取导弹数据和目标数据 dan_list[][3]（编号，x，y） target_list[][4]（编号，x，y，价值）
    // 1.1 读取导弹数据
    const wchar_t* filename = L"MissileParameters8.xls";
    Book* book = xlCreateBook(); // 创建一个 Excel 文档对象
    double dan_list[8][3];
    if (book->load(filename)) { // 加载 Excel 文件
        Sheet* sheet = book->getSheet(0); // 获取第一个工作表
        if (sheet) {
            //int rowCount = sheet->lastRow(); // 获取行数
            //int colCount = sheet->lastCol(); // 获取列数

            for (int i = 0; i < 8; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    dan_list[i][j] = sheet->readNum(i + 1, j); // 读取单元格数据
                }

            }
        }
        else {
            std::cerr << "Error accessing sheet!" << std::endl;
        }
    }
     book->release(); // 释放 Excel 文档对象

     //1.2 读取目标数据——————————————————————————
     const wchar_t* filename2 = L"Target3.xls";
     Book* book2 = xlCreateBook(); // 创建一个 Excel 文档对象
     double target_list[3][4];
     if (book2->load(filename2)) { // 加载 Excel 文件
         Sheet* sheet = book->getSheet(0); // 获取第一个工作表
         if (sheet) {
             //int rowCount = sheet->lastRow(); // 获取行数
             //int colCount = sheet->lastCol(); // 获取列数

             for (int i = 0; i < 3; ++i)
             {
                 for (int j = 0; j < 3; ++j)
                 {
                     target_list[i][j] = sheet->readNum(i + 1, j); // 读取单元格数据
                 }
                 target_list[i][3] = sheet->readNum(i + 1, 8);
             }
         }
         else {
             std::cerr << "Error accessing sheet!" << std::endl;
         }
     }
     book2->release(); // 释放 Excel 文档对象
     //打印list信息
     //int rows = sizeof(target_list) / sizeof(target_list[0]);
     //int cols = sizeof(target_list[0]) / sizeof(target_list[0][0]);
     //printArray(target_list, rows, cols);

     // （2）集成假设数据  list_test 导弹数量，目标数量，弹编号，弹x，弹y，弹角度，，，，，，目标编号，目标x，目标y，目标角度，目标价值
     std::vector<double> list_test;
     int n_dan = 8;
     int n_tar = 1;

     list_test.push_back(n_dan);
     list_test.push_back(n_tar);

     for (int i = 0; i < n_dan; ++i) {
         list_test.push_back(dan_list[i][0]);
         list_test.push_back(dan_list[i][1]);
         list_test.push_back(dan_list[i][2]);
         list_test.push_back(0);
     }

     for (int j = 0; j < n_tar; ++j) {
         list_test.push_back(target_list[j ][0]);
         list_test.push_back(target_list[j][1]);
         list_test.push_back(target_list[j ][2]);
         list_test.push_back(pi);
         list_test.push_back(target_list[j ][3]);
     }

     // 打印list_test
     //std::cout << "Input list:" << std::endl;
     //for (double value : list_test) {
     //    std::cout << value << " ";
     //}
     //Target aaa = Target(0);
     //std::cout << aaa.number << std::endl;

     //Missile missile1(1);
     //std::cout << missile1.number << std::endl;


     std::unordered_map<int, Missile> dan_map;
     std::unordered_map<int, Target> target_map;

     int dan_num, int ship_num = int(list_test[0]), int(list_test[1]);
     int index, bullet_num, double x, y, theta;
     std::unordered_set<int> dan_live;

     for (int i = 0; i < dan_num; ++i)
     {
         index, x, y, theta = list_test[4 * i + 2], list_test[4 * i + 3], list_test[4 * i + 4], list_test[4 * i + 5];
         dan_map[index]=(Missile(index));
         dan_map[index].location[0], dan_map[index].location[1] = x, y;
         dan_map[index].theta = 0;
         dan_live.insert(index);
     }

     for (int i = 0; i < ship_num; ++i) {
         index, x, y, theta , bullet_num = list_test[5 * i+ 2+4* dan_num], list_test[5 * i + 3 + 4 * dan_num], list_test[5 * i + 4 + 4 * dan_num],
                                               list_test[5* i + 5 + 4 * dan_num], list_test[5 * i + 6 + 4 * dan_num];
         target_map[index] = (Target(index));
         target_map[index].location[0], dan_map[index].location[1] = x, y;
         target_map[index].theta = pi;
         target_map[index].value = bullet_num * 100;
         dan_live.insert(index);
     }
     for (int i = 0; i < 2; ++i)
     {
         std::cout << target_map[i].number << std::endl;
     }

    return 0;


}
