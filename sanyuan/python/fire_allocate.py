import copy
import random
import time as code_time
import numpy as np
import math
from math import inf, sin, cos, acos, sqrt, atan2, pi
from data_processing import read_data
from math import pi


def dis(a, b):
    sum2 = 0
    for i in range(len(a)):
        sum2 += (a[i]-b[i])**2
    return sqrt(sum2)


def add2d_dict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})


def del2d_dict(thedict, key_a, key_b):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            del thedict[key_a][key_b]


def rad_normol(theta):
    if theta > pi:
        theta -= 2*pi
    elif theta < -pi:
        theta += 2*pi
    return theta

def dis_lon_lat(lon1, lat1, lon2, lat2):
    # lon经度
    # lat纬度
    R = 6371  # 地球半径（千米）
    # 度转弧度
    lon1rad, lat1rad, lon2rad, lat2rad = lon1*pi/180, lat1*pi/180, lon2*pi/180, lat2*pi/180
    d = R * acos(cos(lat1rad) * cos(lat2rad) * cos(lon1rad - lon2rad) + sin(lat1rad) * sin(lat2rad))
    # 百公里误差数米，千公里误差数十米
    return d


def mod(a, b):
    return a % b


def lsl(alpha, beta, d):
    p_squared = 2 + (d * d) - (2 * cos(alpha - beta)) + (2 * d * (sin(alpha) - sin(beta)))
    if p_squared < 0:
        return [inf, inf, inf, inf]

    tmp0 = d + sin(alpha) - sin(beta)
    tmp1 = atan2((cos(beta) - cos(alpha)), tmp0)
    t = mod((-alpha + tmp1), 2 * pi)
    p = sqrt(p_squared)
    q = mod((beta - tmp1), 2 * pi)
    length = t+p+q

    return [length, t, p, q]


def lsr(alpha, beta, d):
    p_squared = -2 + d * d + 2 * cos(alpha - beta) + 2 * d * (sin(alpha) + sin(beta))
    if p_squared < 0:
        return [inf, inf, inf, inf]
    p = sqrt(p_squared)
    tmp2 = atan2((-cos(alpha) - cos(beta)), (d + sin(alpha) + sin(beta))) - atan2(-2.0, p)
    t = mod((-alpha + tmp2), 2 * pi)
    q = mod((-mod(beta, 2 * pi) + tmp2), 2 * pi)
    length = t + p + q

    return [length, t, p, q]


def rsl(alpha, beta, d):
    p_squared = (d * d) - 2 + (2 * cos(alpha - beta)) - (2 * d * (sin(alpha) + sin(beta)))
    if p_squared < 0:
        return [inf, inf, inf, inf]
    p = sqrt(p_squared)
    tmp2 = atan2((cos(alpha) + cos(beta)), (d - sin(alpha) - sin(beta))) - atan2(2.0, p)
    t = mod((alpha - tmp2), 2 * pi)
    q = mod((beta - tmp2), 2 * pi)
    length = t + p + q

    return [length, t, p, q]


def rsr(alpha, beta, d):
    p_squared = 2 + (d * d) - (2 * cos(alpha - beta)) + (2 * d * (sin(beta) - sin(alpha)))
    if p_squared < 0:
        return [inf, inf, inf, inf]
    tmp0 = d - sin(alpha) + sin(beta)
    tmp1 = atan2((cos(alpha) - cos(beta)), tmp0)
    t = mod((alpha - tmp1), 2 * pi)
    p = sqrt(p_squared)
    q = mod((-beta + tmp1), 2 * pi)
    length = t + p + q

    return [length, t, p, q]


def rlr(alpha, beta, d):
    tmp_rlr = (6. - d * d + 2 * cos(alpha - beta) + 2 * d * (sin(alpha) - sin(beta))) / 8.
    if abs(tmp_rlr) > 1:
        return [inf, inf, inf, inf]
    p = mod((2 * pi - acos(tmp_rlr)), 2 * pi)
    t = mod((alpha - atan2(cos(alpha) - cos(beta), d - sin(alpha) + sin(beta)) + mod(p / 2, 2 * pi)), 2 * pi)
    q = mod((alpha - beta - t + mod(p, 2 * pi)), 2 * pi)
    length = t + p + q

    return [length, t, p, q]


def lrl(alpha, beta, d):
    tmp_lrl = (6. - d * d + 2 * cos(alpha - beta) + 2 * d * (- sin(alpha) + sin(beta))) / 8.
    if abs(tmp_lrl) > 1:
        return [inf, inf, inf, inf]
    p = mod((2 * pi - acos(tmp_lrl)), 2 * pi)
    t = mod((-alpha - atan2(cos(alpha) - cos(beta), d + sin(alpha) - sin(beta)) + p / 2), 2 * pi)
    q = mod((mod(beta, 2 * pi) - alpha - t + mod(p, 2 * pi)), 2 * pi)
    length = t + p + q

    return [length, t, p, q]


def dubins_min_len(p1, theta1, p2, theta2, r):
    # 经纬度距离转化，建立相对平面直角坐标系
    # dx = dis_lon_lat(p1[0], p2[1], p2[0], p2[1])
    # dy = dis_lon_lat(p2[0], p1[1], p2[0], p2[1])

    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    d = sqrt(dx ** 2 + dy ** 2) / r

    theta = mod(atan2(dy, dx), 2 * pi)
    alpha = mod((theta1 - theta), 2 * pi)
    beta = mod((theta2 - theta), 2 * pi)

    length = lsl(alpha, beta, d)
    tmp = lsr(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp
    tmp = rsl(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp

    tmp = rsr(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp

    tmp = rlr(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp

    tmp = lrl(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp

    return r*length[0]


def dubins_segment(seg_param, seg_init, seg_type):
    seg_end = [0, 0, 0]
    if seg_type == 'L':
        seg_end[0] = seg_init[0] + sin(seg_init[2] + seg_param) - sin(seg_init[2])
        seg_end[1] = seg_init[1] - cos(seg_init[2] + seg_param) + cos(seg_init[2])
        seg_end[2] = seg_init[2] + seg_param
    elif seg_type == 'R':
        seg_end[0] = seg_init[0] - sin(seg_init[2] - seg_param) + sin(seg_init[2])
        seg_end[1] = seg_init[1] + cos(seg_init[2] - seg_param) - cos(seg_init[2])
        seg_end[2] = seg_init[2] - seg_param
    elif seg_type == 'S':
        seg_end[0] = seg_init[0] + cos(seg_init[2]) * seg_param
        seg_end[1] = seg_init[1] + sin(seg_init[2]) * seg_param
        seg_end[2] = seg_init[2]
    return seg_end


def dubins_end_point(p1, theta1, p2, theta2, r, s):
    # 经纬度距离转化
    # dx = dis_lon_lat(p1[0], p2[1], p2[0], p2[1])
    # dy = dis_lon_lat(p2[0], p1[1], p2[0], p2[1])
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    d = sqrt(dx ** 2 + dy ** 2) / r

    theta = mod(atan2(dy, dx), 2 * pi)
    alpha = mod((theta1 - theta), 2 * pi)
    beta = mod((theta2 - theta), 2 * pi)

    length = lsl(alpha, beta, d)
    type_d = 'LSL'
    tmp = lsr(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp
        type_d = 'LSR'
    tmp = rsl(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp
        type_d = 'RSL'
    tmp = rsr(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp
        type_d = 'RSR'
    tmp = rlr(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp
        type_d = 'RLR'
    tmp = lrl(alpha, beta, d)
    if tmp[0] < length[0]:
        length = tmp
        type_d = 'LRL'

    seg_param = s / r
    seg_init = [0, 0, theta1]
    if seg_param <= length[1]:
        seg_end = dubins_segment(seg_param, seg_init, type_d[0])
        x = p1[0] + seg_end[0] * r
        y = p1[1] + seg_end[1] * r
        theta = mod(seg_end[2], 2*pi)
        return x, y, theta

    mid1 = dubins_segment(length[1], seg_init, type_d[0])
    if seg_param <= length[1]+length[2]:
        seg_end = dubins_segment(seg_param - length[1], mid1, type_d[1])
        x = p1[0] + seg_end[0] * r
        y = p1[1] + seg_end[1] * r
        theta = mod(seg_end[2], 2 * pi)
        return x, y, theta

    mid2 = dubins_segment(length[2], mid1, type_d[1])
    seg_end = dubins_segment(seg_param - length[1]-length[2], mid2, type_d[2])
    x = p1[0] + seg_end[0] * r
    y = p1[1] + seg_end[1] * r
    theta = mod(seg_end[2], 2 * pi)
    return x, y, theta

#############################################################env

class Message:
    def __init__(self, receiver, sender, time, text):
        self.receiver = receiver  # 接收方编号
        self.sender = sender  # 发送方编号
        self.send_time = time  # 当前时间
        self.receive_time = 0  # 在环境中传播后，接收到的时间

        # 环境时延
        # delay = 3 + 3 * sin(0.05 * time + 1) + 2 * sin(1.22 * time + .5 * time + 0.2) + 1.5 * sin(2.2 * time + 1) + cos(
        #     1.8 * time + 1) ** 2
        # if delay < 1:
        #     delay = 1
        delay = 1
        self.receive_time = time + delay
        self.delay = delay

        # 信息内容
        self.text = text


class Missile:
    def __init__(self, number):
        self.number = number
        self.location = np.array([0., 0., 0.])
        self.s_range = 400.0  # 剩余航程 km
        self.velocity = 0.200  # 速度 km/s
        self.radius = 20  # 转弯半径 km
        self.theta = 0  # 航向角 rad
        self.hit_theta = 0  # 打击入射角 rad
        self.value = 100  # 价值
        self.commu_distance = 20  # 通信距离 km
        self.commu_prob = 1  # 通信质量 0~1之间
        self.health = 1  # 当前状态 0 损毁， 1 健康
        self.hit = 0  # 是否击中目标， 0 未击中， 1 击中
        self.angle_min = pi / 60  # 相邻弹最小打击角度间隔，3度
        self.task1 = Target(0)  # 执行任务信息
        self.task2 = Target(0)  # 携带任务信息
        self.advantage_val = 0  # 打击优势值
        self.hit_distance = 0  # 距离目标剩余打击距离
        self.hit_theta_dict = {}  # 用于存放打击同一任务成员的打击角度值,key 1 为任务编号,key 2为导弹编号,值为打击角度
        self.hit_theta_map = []  # 用于存放打击同一任务成员的打击角度值
        self.hit_theta_left = [-10, 0, 0]  # 用于存放打击同一任务成员的打击角度值紧临的小于值,存储格式[角度，导弹编号，任务编号]
        self.hit_theta_right = [10, 0, 0]  # 用于存放打击同一任务成员的打击角度值紧临的大于值,存储格式[角度，导弹编号，任务编号]
        self.task_value_dict = {}  # 用于存放打击任务成员的剩余航程优势函数值，key 1 为任务编号,key 2为导弹编号
        self.task_value = 0  # 执行任务的优势值，剩余航程

        self.phase = 0  # 通信阶段 0 为第一次通信，1为发送决策内容，2为决策确认
        self.message = []  # 报文内容

        self.connect_buffer = {}  # 连接缓存器，存放保持连接的飞行器，用飞行器编号作为索引
        self.state_time = {}  # 最新接受到飞行器消息的时间，用飞行器编号作为索引，用于判断飞行器是否存活
        self.connect_time = {}  # 连接时间，飞行器状态最新更新时间，用飞行器编号作为索引
        self.task_buffer = []  # 从盟友回收的任务，或侦察获得的额外任务
        self.interrupt_buffer = {}  # 中断缓存器
        self.hit_buffer = {}  # 打击角缓存器，用于存在打击相同目标弹群的打击角，键为导弹编号，值为打击角度；按照打击角度进行排序，便于解除冲突

        self.ally_buffer = {}  # 协作联盟，用于存放盟友的编号，盟友的信息数据在connect_buffer中，键为盟友编号，值为盟友状态，0为正常，
        # 1 2 3为发送脱离信息次数
        self.ally_num = 0  # 盟友个数
        self.ally_tmp_number = 0  # 待确认阶段
        self.ally_stage = 0  # 1：请求加入，2：回复请求，3：确认执行
        self.ally_time = 0  # 计时器
        self.ally_flag = 0  # 请求标识
        self.ally_deadline = 10  # 通信等待截止时间
        self.ally_reply_buffer = set()  # 等待回复
        self.ally_reply_timer = 0  # 等待回复计时器
        self.ally_reply_flag = 0  # 等待计时标志
        self.ally_reply_message = 0  # 需要重发的信息
        self.interrupt_time = 10  # 失联时间，如果超过interrupt_time，则主动传呼盟友；超过2*interrupt_time，则回收盟友任务

        # self.request_buffer = {}  # 请求缓存器，key 为飞行器编号，内容为[flag,time,flag,time],【发出的请求，发出时间，接收到的请求，接收时间】
        # self.receive_message_buffer = {}  # 接收到的内容

        self.message_buffer = []  # 用于暂放消息

        self.buy_stage = 0  # 1:广播任务，2：投标，3：中标结果，4：确认执行
        self.buy_flag = 0  # 标识初始化
        self.buy_time = 0  # 计时器初始化
        self.buy_number = 0  # 存放交易对象编号
        self.buy_need = 0  # 买入需求:1为有需求，0为无需求

        self.buy_time_begin = 2
        self.buy_deadline = 10
        self.buy_target = Target(0)
        self.buy_bidder_buffer = {}  # 存放标书

        self.buy_task_buffer = {}  # 存放接收到的任务，安照发送方编号存储
        self.buy_task_value = {}  # 存放接收到的任务优势值，安照发送方编号存储
        self.buy_reply_buffer = set()  # 等待回复
        self.buy_reply_timer = 0  # 等待回复计时器
        self.buy_reply_flag = 0  # 等待计时标志
        self.buy_reply_message = 0  # 需要重发的信息
        # self.buy_connect_time = {}  # 连接时间，用飞行器编号作为索引

        self.sell_stage = 0
        self.sell_time = 0
        self.sell_flag = 0  # 0:买，1：卖
        self.sell_need = 0
        self.sell_num = 0  # 交易对象编号
        self.sell_bullet_num = 0
        self.sell_reply_buffer = set()  # 等待回复
        self.sell_reply_timer = 0  # 等待回复计时器
        self.sell_reply_flag = 0  # 等待计时标志
        self.sell_reply_message = 0  # 需要重发的信息

        self.exchange_number = 0
        self.exchange_stage = 0  # 1：请求交换，2：回复请求，3：确认执行
        self.exchange_time = 0
        self.exchange_flag = 0
        self.exchange_deadline = 10
        self.exchange_target = Target(0)
        self.exchange_reply_buffer = set()  # 等待回复
        self.exchange_reply_timer = 0  # 等待回复计时器
        self.exchange_reply_flag = 0  # 等待计时标志
        self.exchange_reply_message = 0  # 需要重发的信息
        self.target_buffer = []  # 暂存用于交换的目标
        self.exchange_text = []  # 用于存放交换的信息

    def task_reset(self):
        self.hit_theta_left = [-10, 0, 0]  # 用于存放打击同一任务成员的打击角度值紧临的小于值,存储格式[角度，导弹编号，任务编号]
        self.hit_theta_right = [10, 0, 0]  # 用于存放打击同一任务成员的打击角度值紧临的大于值,存储格式[角度，导弹编号，任务编号]

    def buy_reset(self):
        self.buy_stage = 0  # 1:广播任务，2：投标，3：中标结果，4：确认执行
        self.buy_flag = 0  # 标识初始化
        self.buy_time = 0  # 计时器初始化
        self.buy_number = 0  # 存放交易对象编号
        self.buy_need = 0  # 买入需求:1为有需求，0为无需求
        self.buy_target.number = 0  # 交易任务初始化

    def sell_reset(self):
        self.sell_stage = 0
        self.sell_time = 0
        self.sell_flag = 0  # 0:买，1：卖

    def ally_reset(self):
        self.ally_tmp_number = 0  # 待确认阶段
        self.ally_stage = 0  # 1：请求加入，2：回复请求，3：确认执行
        self.ally_time = 0  # 计时器
        self.ally_flag = 0  # 请求标识

    def exchange_reset(self):
        self.exchange_number = 0
        self.exchange_stage = 0  # 1：请求交换，2：回复请求，3：确认执行
        self.exchange_time = 0
        self.exchange_flag = 0

    def update_missile(self, location, s_range, commu_prob, health, task1, task2):
        self.location = location
        self.s_range = s_range  # 航程 km
        self.commu_prob = commu_prob  # 通信质量 0~1之间
        self.health = health  # 当前状态 0 损毁， 1 健康
        self.task1 = task1  # 执行任务编号
        self.task2 = task2  # 携带任务编号

    def advantage(self, ship, hit_theta=None, distance=None):
        print('原始的优势值')
        # 距离优势
        if hit_theta is None and distance is None:
            hit_theta = ship.theta_domain0[0] + (ship.theta_domain0[1] - ship.theta_domain0[0]) * random.random()
            hit_theta_abs = rad_normol(ship.theta + hit_theta)  # 绝对角度
            distance = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs, self.radius)

            hit_theta_tmp = ship.theta_domain1[0] + (ship.theta_domain1[1] - ship.theta_domain1[0]) * random.random()
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp)
            distance_tmp = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs_tmp, self.radius)
            if distance_tmp < distance:
                hit_theta = hit_theta_tmp
                distance = distance_tmp
            hit_theta_tmp = ship.theta_domain2[0] + (ship.theta_domain2[1] - ship.theta_domain2[0]) * random.random()
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp)
            distance_tmp = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs_tmp, self.radius)
            if distance_tmp < distance:
                hit_theta = hit_theta_tmp
                distance = distance_tmp
            hit_theta_tmp = ship.theta_domain3[0] + (ship.theta_domain3[1] - ship.theta_domain3[0]) * random.random()
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp)
            distance_tmp = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs_tmp, self.radius)
            if distance_tmp < distance:
                hit_theta = hit_theta_tmp
                distance = distance_tmp
        elif distance is None:
            hit_theta_abs = rad_normol(ship.theta + hit_theta)  # 绝对角度
            distance = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs, self.radius)
        # d_bar = self.location - ship.location
        # distance = sqrt(d_bar[0] ** 2 + d_bar[1] ** 2)
        if distance > self.s_range:
            return 0
        elif self.s_range < 50:
            vd = 1 + 500 / distance
        else:
            vd = 1 + distance / self.s_range + self.s_range / distance

        # 代价优势
        vv = ship.value / self.value

        return vd * vv

    def decision_making(self, ship, hit_theta=None):
        # 距离优势
        if hit_theta is None:
            hit_theta = ship.theta_domain0[0] + (ship.theta_domain0[1] - ship.theta_domain0[0]) * random.random()
            hit_theta_abs = rad_normol(ship.theta + hit_theta)  # 绝对角度
            distance = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs, self.radius)

            hit_theta_tmp = ship.theta_domain1[0] + (ship.theta_domain1[1] - ship.theta_domain1[0]) * random.random()
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp)
            distance_tmp = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs_tmp, self.radius)
            if distance_tmp < distance:
                hit_theta = hit_theta_tmp
                distance = distance_tmp
            hit_theta_tmp = ship.theta_domain2[0] + (ship.theta_domain2[1] - ship.theta_domain2[0]) * random.random()
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp)
            distance_tmp = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs_tmp, self.radius)
            if distance_tmp < distance:
                hit_theta = hit_theta_tmp
                distance = distance_tmp
            hit_theta_tmp = ship.theta_domain3[0] + (ship.theta_domain3[1] - ship.theta_domain3[0]) * random.random()
            hit_theta_abs_tmp = rad_normol(ship.theta + hit_theta_tmp)
            distance_tmp = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs_tmp, self.radius)
            if distance_tmp < distance:
                hit_theta = hit_theta_tmp
                distance = distance_tmp
        else:
            hit_theta_abs = rad_normol(ship.theta + hit_theta)  # 绝对角度
            distance = dubins_min_len(self.location, self.theta, ship.location, hit_theta_abs, self.radius)
        # d_bar = self.location - ship.location
        # distance = sqrt(d_bar[0] ** 2 + d_bar[1] ** 2)
        if distance > self.s_range:
            return 0, 0, 0
        elif self.s_range < 50:
            vd = 1 + 500 / distance
        else:
            vd = 1 + distance / self.s_range + self.s_range / distance

        # 代价优势
        vv = ship.value / self.value

        return vd * vv, hit_theta, distance

    def moving1(self, dt):

        s = self.velocity * dt
        x, y, theta = dubins_end_point(self.location, self.theta, self.task1.location,
                                       self.task1.theta + self.hit_theta,
                                       self.radius, s)
        self.location[0] = x
        self.location[1] = y
        self.theta = theta

    def moving2(self, dt):
        if self.task1.number > 0:  # 有打击目标
            val, hit_theta, distance = self.decision_making(self.task1)
            self.update_hit_theta_map()
            if distance < self.hit_distance and self.check_new_hit_theta(hit_theta):
                self.hit_theta = hit_theta
                self.hit_distance = distance
                # self.task_reset()
            self.angle_check()  # 检查是否存在角度冲突

            # s = self.velocity * dt
            # self.hit_distance -= s
            self.hit_distance = distance
            self.advantage_val = self.advantage(self.task1, distance=self.hit_distance)

    def moving(self, dt):  # 飞行器运动
        if self.task1.number > 0:  # 有打击目标
            val, hit_theta, distance = self.decision_making(self.task1)
            self.update_hit_theta_map()
            if distance < self.hit_distance and self.check_new_hit_theta(hit_theta):
                self.hit_theta = hit_theta
                self.hit_distance = distance
                # self.task_reset()
            self.angle_check()  # 检查是否存在角度冲突

        s = self.velocity * dt
        self.hit_distance -= s
        self.advantage_val = self.advantage(self.task1, distance=self.hit_distance)

        x, y, theta = dubins_end_point(self.location, self.theta, self.task1.location,
                                       self.task1.theta + self.hit_theta,
                                       self.radius, s)
        self.location[0] = x
        self.location[1] = y
        self.theta = theta
        d = self.task1.location - self.location
        s3 = sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
        s2 = sqrt(d[0] ** 2 + d[1] ** 2)
        # if s3 < 2:
        #     self.health = 0
        # self.hit = 1
        # print('击中： 目标', self.task1.number)
        if s2 < 10 and ((self.task1.number != 0 and self.task2.number == 0) or self.radius * pi * 2 > self.s_range):
            # 下降高度，攻击目标
            th = s2 / self.velocity
            self.location[2] -= self.location[2] / th * dt

        self.s_range -= self.velocity * dt
        s2 = dubins_min_len(self.location, self.theta, self.task1.location, 0, self.radius)
        self.task_value = self.s_range - s2  # 剩余飞行里程
        if self.s_range < 0:
            self.health = 0
            # print('飞行器 %d 燃料耗尽' % self.number)

    def ageing(self, time):
        # 时效性管理
        for num in list(self.connect_time.keys()):
            if time - self.connect_time[num] > 10:
                if self.task1.number and self.task1.number == self.connect_buffer[num].task1.number:
                    del2d_dict(self.task_value_dict, self.task1.number, self.connect_buffer[num].number)
                self.connect_time.pop(num)
        # 执行相同任务成员时效性管理
        if self.task1.number and self.task1.number in self.task_value_dict:
            for num in list(self.task_value_dict[self.task1.number]):
                if time - self.task_value_dict[self.task1.number][num][1] > 5 or \
                        self.task_value_dict[self.task1.number][num][0] < 0:
                    del2d_dict(self.task_value_dict, self.task1.number, num)
        #         if num in self.request_buffer:
        #             self.request_buffer.pop(num)
        #         if num in self.receive_message_buffer:
        #             self.receive_message_buffer.pop(num)
        if self.buy_flag:
            self.buyer_age()
        if self.exchange_flag:
            self.exchange_age()
        if self.sell_flag:
            self.seller_age()
        if self.ally_flag:
            self.ally_age()
        if self.task1.number == 0 or self.task2.number == 0:
            self.buy_need = 1
        else:
            self.buy_need = 0
        if self.task2.number == 0:
            self.sell_need = 0
            if self.target_buffer:
                copy_target(self.task2, self.target_buffer.pop())
        else:
            self.sell_need = 1

        # 目标重复自查
        self.call_off()
        # 目标回收
        if self.target_buffer and self.task1.number == 0:
            val = 0
            num = -1
            for i in range(len(self.task_buffer)):
                tmp_val, tmp_hit_theta, tmp_distance = self.decision_making(self.target_buffer[-1])
                if tmp_val > val:
                    val = tmp_val
                    num = i
                    hit_theta = tmp_hit_theta
                    distance = tmp_distance
                    break
            if num > -1:
                self.task1 = self.target_buffer.pop(num)
                self.advantage_val = val
                self.hit_theta = hit_theta
                self.hit_distance = distance
                self.task_reset()
        if self.target_buffer and self.task2.number == 0:  # 将回收任务进行传播
            self.task2 = self.target_buffer.pop()

    # 重传函数
    def reply(self):
        message = []
        if self.exchange_reply_flag:
            self.exchange_reply_timer += 1
            if self.exchange_reply_timer == 2:
                message.append(self.exchange_reply_message)
            elif self.exchange_reply_timer == 5:
                message.append(self.exchange_reply_message)
                self.exchange_reply_message = 0
                self.exchange_reply_timer = 0
                self.exchange_reply_flag = 0

        if self.buy_reply_flag:
            self.buy_reply_timer += 1
            if self.buy_reply_timer == 2:
                message.append(self.buy_reply_message)
            elif self.buy_reply_timer == 5:
                message.append(self.buy_reply_message)
                self.buy_reply_message = 0
                self.buy_reply_timer = 0
                self.buy_reply_flag = 0
        if self.ally_reply_flag:
            self.ally_reply_timer += 1
            if self.ally_reply_timer == 2:
                message.append(self.ally_reply_message)
            elif self.ally_reply_timer == 5:
                message.append(self.ally_reply_message)
                self.ally_reply_message = 0
                self.ally_reply_timer = 0
                self.ally_reply_flag = 0
        # if self.message_buffer:  # 将缓存器中数据发送出去
        #     for mess in self.message_buffer:
        #         message.append(mess)
        for _ in range(len(self.message_buffer)):
            message.append(self.message_buffer.pop())
        return message

    def buyer_age(self):  # 买方时效性评估
        self.buy_time += 1
        if self.buy_time > self.buy_deadline:
            self.buy_flag = 0
            self.buy_stage = 0
            self.buy_time = 0
            self.buy_number = 0

    def seller_age(self):
        self.sell_time += 1
        if self.sell_time > self.buy_deadline:
            self.sell_flag = 0
            self.sell_stage = 0
            self.sell_time = 0

    def exchange_age(self):
        self.exchange_time += 1
        if self.exchange_time > self.exchange_deadline:
            self.exchange_time = 0
            self.exchange_stage = 0
            self.exchange_flag = 0

    def ally_age(self):
        self.ally_time += 1
        if self.ally_time > self.ally_deadline:
            self.ally_time = 0
            self.ally_stage = 0
            self.ally_flag = 0

    # 买卖部分：
    def broadcast(self, dan2_number, time):
        if self.sell_stage == 0 and self.task2.number > 0:  # 卖的需求
            self.sell_bullet_num = self.task2.bullet_num  # 0变1
            self.sell_num = 0  # 交易对象编号
            self.sell_flag = 1  # 0:买，1：卖

            # 广播任务
            text = [1, 1, self.task2]  # 任务类型，任务阶段，任务信息

            return Message(dan2_number, self.number, time, text)  # 生成消息
        return 0

    def buy_doing(self, time):
        # 根据招标信息，投标
        if self.buy_flag == 0:
            self.buy_time += 1
            if self.buy_time >= self.buy_time_begin:
                self.buy_time = 0
                if self.buy_task_value:
                    # 任务评价
                    a1 = sorted(self.buy_task_value.items(), key=lambda x: x[1])  # 按照优势值排序
                    self.buy_number = a1[0][0]
                    copy_target(self.buy_target, self.buy_task_buffer[self.buy_number])
                    task1_flag, task2_flag = 0, 0
                    val, hit_theta, hit_distance = self.buy_task_value[self.buy_number]
                    if self.task1.number == 0 and val > 0:
                        task1_flag = self.buy_target.number
                        self.advantage_val = val
                        self.hit_theta = hit_theta
                        self.hit_distance = hit_distance
                        self.buy_target.bullet_num -= 1
                    if self.task2.number == 0 and self.buy_target.bullet_num > 1:
                        task2_flag = self.buy_target.number

                    # 清空
                    self.buy_task_value = {}
                    self.buy_task_buffer = {}

                    if task1_flag or task2_flag:
                        self.buy_reply_buffer.add(self.buy_number)
                        self.buy_stage = 2
                        self.buy_flag = 1
                        text = [1, 2, task1_flag, task2_flag]
                        return Message(self.buy_number, self.number, time, text)
                    else:
                        return 0
        return 0

    def bidder(self, time):
        # 任务传播结果确认（任务执行的买卖已经执行）
        self.sell_time += 1
        if self.sell_time > self.buy_deadline:
            self.sell_stage = 3
            self.sell_time = 0

        if self.sell_stage == 3:
            if self.buy_bidder_buffer:
                mess = []
                for num in list(self.buy_bidder_buffer.keys()):
                    message = self.buy_bidder_buffer.pop(num)
                    task1_flag, task2_flag = 0, 0
                    if message.text[2] == self.task2.number and self.sell_bullet_num > 0:
                        task1_flag = self.task2.number
                        self.sell_bullet_num -= 1
                    if message.text[3] == self.task2.number and self.sell_bullet_num > 0:
                        task2_flag = self.task2.number
                        self.sell_bullet_num -= 1
                    text = [1, 3, task1_flag, task2_flag]  # 1, 阶段3， task1成交任务编号，task2成交任务编号
                    mess.append(Message(message.sender, self.number, time, text))
                    if task1_flag or task2_flag:
                        self.sell_num += 1
                        self.sell_reply_buffer.add(message.sender)  # 等待此对象回复
                return mess
        return 0

    # def doing(self):  # 执行决定
    #     for message in self.message_buffer:
    #         pass

    # 联盟部分：
    def ally_request(self, time):
        if self.ally_num == 0 and self.ally_stage == 0 and self.task1.number != 0:
            # 需要加入联盟，从建立连接的成员中选择合适的建立联盟
            best_num = 0
            val = 1e8
            for num in list(self.connect_buffer.keys()):
                if self.connect_buffer[num].ally_num < 2:
                    if self.connect_buffer[num].task1.number == 0:
                        temp = 10
                    else:
                        temp = dis(self.task1.location, self.connect_buffer[num].task1.location)
                    if temp < val:
                        best_num = num
                        val = temp
            if best_num != 0:
                # 加入联盟请求
                text = [3, 1, 1]  # 任务类型3，任务阶段1，任务信息:0 拒绝，1同意
                self.ally_stage = 1
                self.ally_flag = 1
                self.ally_tmp_number = best_num
                self.message_buffer.append(Message(best_num, self.number, time, text))  # 生成消息

        if self.ally_num > 0:
            for num in list(self.ally_buffer.keys()):
                if num not in self.connect_buffer:  # 不在连接器中，询问其具体信息
                    text = [3, 0, 1]
                    self.message_buffer.append(Message(num, self.number, time, text))
                else:
                    if dis(self.location, self.connect_buffer[num].location) > 0.8 * self.commu_distance:
                        if self.ally_buffer[num] < 4:
                            self.ally_buffer[num] += 1
                        if self.ally_buffer[num] == 3:
                            self.ally_buffer.pop(num)  # 脱离联盟
                            self.ally_num -= 1
                        text = [3, 1, 0]  # 连续发送三次发送脱离信息
                        self.message_buffer.append(Message(num, self.number, time, text))
                    else:
                        self.ally_buffer[num] = 0

    def ally(self, time):
        # 盟友状态检查
        for num in list(self.ally_buffer.keys()):
            if time - self.state_time[num] > 5:
                self.message_buffer.append(Message(num, self.number, time, [3, 0, 1]))  # 发送信息询问状态
            if time - self.state_time[num] > 10:  # 判定为失联,回收任务
                if num not in self.connect_buffer:  # 不在连接器中，删除其具体信息
                    self.ally_buffer.pop(num)  # 移除联盟
                    self.ally_num -= 1  # 数量减一
                else:
                    dan = self.connect_buffer[num]
                    if dan.task1.number != 0:
                        self.target_buffer.append(dan.task1)  # 判定为失联,回收任务
                    if dan.task2.number != 0:
                        self.target_buffer.append(dan.task2)  # 判定为失联,回收任务
                    self.interrupt_buffer[dan.number] = dan  # 添加到失联缓存器中
                    self.ally_buffer.pop(num)  # 移除联盟
                    self.ally_num -= 1  # 数量减一

    def call_off(self):
        # if self.task1.number and self.task_value < 0:
        #     self.task1.number = 0  # 取消任务执行
        if self.task1.number and self.task1.number in self.task_value_dict and self.task1.bullet_real <= \
                len(self.task_value_dict[self.task1.number]):
            # 自查是否取消任务
            flag = 0
            order = 1
            for num in list(self.task_value_dict[self.task1.number]):
                if self.task_value_dict[self.task1.number][num][0] < 0:
                    del2d_dict(self.task_value_dict, self.task1.number, num)

                if self.task_value_dict[self.task1.number][num][0] < self.task_value:
                    order += 1
                    if order > self.task1.bullet_real:
                        flag = 1
                        break
            # if flag:
            #     self.task1.number = 0  # 取消任务执行

    def angle_check(self):
        if self.hit_theta_left[0] > self.hit_theta or self.hit_theta_right[0] < self.hit_theta:
            self.task_reset()
        delta1 = self.hit_theta - self.hit_theta_left[0]
        delta2 = self.hit_theta_right[0] - self.hit_theta
        if delta1 < self.angle_min or delta2 < self.angle_min:  # 存在角度过小，需要角度自调
            if self.hit_theta_right[0] - self.hit_theta_left[0] < 2 * self.angle_min:  # 两边都紧，往中间调整
                self.hit_theta = (self.hit_theta_left[0] + self.hit_theta_right[0]) / 2
            elif delta1 < self.angle_min:  # 左边紧，向右调整
                self.hit_theta = self.hit_theta_left[0] + self.angle_min
            else:  # 右边紧，向左调整
                self.hit_theta = self.hit_theta_right[0] - self.angle_min
            self.angle_domain()  # 可行域检查

    def angle_domain(self):
        self.hit_theta = rad_normol(self.hit_theta)  # 正则化到【-pi，pi】
        if self.hit_theta > 0:  # 左侧区
            if self.hit_theta > 1.57:  # 左后侧区
                if self.hit_theta < self.task1.theta_domain3[0]:
                    self.hit_theta = self.task1.theta_domain3[0] + 0.02  # 调整到可行区，并1度余量
                elif self.hit_theta > self.task1.theta_domain3[1]:
                    self.hit_theta = self.task1.theta_domain3[1] - 0.02  # 调整到可行区，并1度余量
            else:  # 左前侧区
                if self.hit_theta < self.task1.theta_domain2[0]:
                    self.hit_theta = self.task1.theta_domain2[0] + 0.02  # 调整到可行区，并1度余量
                elif self.hit_theta > self.task1.theta_domain2[1]:
                    self.hit_theta = self.task1.theta_domain2[1] - 0.02  # 调整到可行区，并1度余量
        else:  # 右侧区
            if self.hit_theta < -1.57:  # 右后侧区
                if self.hit_theta < self.task1.theta_domain0[0]:
                    self.hit_theta = self.task1.theta_domain0[0] + 0.02  # 调整到可行区，并1度余量
                elif self.hit_theta > self.task1.theta_domain0[1]:
                    self.hit_theta = self.task1.theta_domain0[1] - 0.02  # 调整到可行区，并1度余量
            else:  # 右前侧区
                if self.hit_theta < self.task1.theta_domain1[0]:
                    self.hit_theta = self.task1.theta_domain1[0] + 0.02  # 调整到可行区，并1度余量
                elif self.hit_theta > self.task1.theta_domain1[1]:
                    self.hit_theta = self.task1.theta_domain1[1] - 0.02  # 调整到可行区，并1度余量

    def update_hit_theta_map(self):
        if self.task1.number in self.hit_theta_dict:
            self.hit_theta_map = sorted(self.hit_theta_dict[self.task1.number].items(), key=lambda x: x[1])
            for i in range(len(self.hit_theta_map)):
                if self.hit_theta_map[i][1][0] > self.hit_theta:
                    self.hit_theta_right = [self.hit_theta_map[i][1][0], self.hit_theta_map[i][0], self.task1.number]
                    if i > 0:
                        self.hit_theta_left = [self.hit_theta_map[i - 1][1][0], self.hit_theta_map[i - 1][0],
                                               self.task1.number]
                    else:
                        self.hit_theta_left = [-10, 0, 0]
                    break
                if i == len(self.hit_theta_map):
                    self.hit_theta_left = [self.hit_theta_map[i][1][0], self.hit_theta_map[i][0], self.task1.number]
                    self.hit_theta_right = [10, 0, 0]

    def check_new_hit_theta(self, new_hit_theta):
        hit_theta_left = [-10, 0, 0]
        hit_theta_right = [10, 0, 0]
        for i in range(len(self.hit_theta_map)):
            if self.hit_theta_map[i][1][0] > new_hit_theta:
                hit_theta_right = [self.hit_theta_map[i][1][0], self.hit_theta_map[i][0], self.task1.number]
                if i > 0:
                    hit_theta_left = [self.hit_theta_map[i - 1][1][0], self.hit_theta_map[i - 1][0], self.task1.number]
                break
            if i == len(self.hit_theta_map) - 1:
                hit_theta_left = [self.hit_theta_map[i][1][0], self.hit_theta_map[i][0], self.task1.number]

        if hit_theta_right[0] - hit_theta_left[0] >= 2 * self.angle_min:
            self.hit_theta_left = hit_theta_left
            self.hit_theta_right = hit_theta_right
            return True
        else:
            return False


class Target:
    def __init__(self, number):
        self.number = number
        self.location = np.array([80., 20., 0.])
        self.theta = pi  # 航向角（弧度）
        self.velocity = 0.016  # 速度 km/s (约60km/h)
        self.value = 0  # 价值
        self.health = 1  # 当前状态 0 损毁， 1 健康
        self.bullet_num = 0  # 任务分配剩余最小饱和打击量
        self.bullet_real = 0  # 饱和打击剩余打击量

        self.theta_head = pi / 6  # 船正前方电子对抗防御区单侧最大夹角
        self.theta_side = pi / 6  # 船侧舷密集火炮防御区单侧最大夹角
        self.theta_end = pi / 9  # 船尾单侧最大夹角
        # 4个可打击角度区间，为相对船头的角度，逆时针为正
        self.theta_domain0 = [- pi + self.theta_end, -pi / 2 - self.theta_side]  # 右后侧区
        self.theta_domain1 = [-pi / 2 + self.theta_side, -self.theta_head]  # 右前侧区
        self.theta_domain2 = [self.theta_head, pi / 2 - self.theta_side]  # 左前侧区
        self.theta_domain3 = [pi / 2 + self.theta_side, pi - self.theta_end]  # 左后侧区

    def update_target(self, location, bullet_num, health):
        self.location = location
        self.health = health  # 当前状态 0 损毁， 1 健康
        self.bullet_num = bullet_num  # 最小饱和打击量

    def moving(self, dt):
        self.location[0] -= self.velocity * dt  # 运动

##########################################tools#######################################################################
def copy_target(new_target: Target, target: Target):
    new_target.number = target.number
    new_target.location[0] = target.location[0]
    new_target.location[1] = target.location[1]
    new_target.value = target.value  # 价值
    new_target.health = target.health  # 当前状态 0 损毁， 1 健康
    new_target.bullet_num = target.bullet_num  # 最小饱和打击量
    new_target.bullet_real = target.bullet_real  # 真实剩余最小饱和打击量


def copy_missile(missile: Missile):
    new_missile = Missile(missile.number)
    new_missile.location[0:3] = missile.location[0:3]
    new_missile.s_range = missile.s_range
    new_missile.commu_prob = missile.commu_prob
    new_missile.task1.number = missile.task1.number
    new_missile.task1.location[0:3] = missile.task1.location[0:3]
    new_missile.task1.value = missile.task1.value
    new_missile.task1.bullet_num = missile.task1.bullet_num
    new_missile.task2.number = missile.task2.number
    new_missile.task2.location[0:3] = missile.task2.location[0:3]
    new_missile.task2.value = missile.task2.value
    new_missile.task2.bullet_num = missile.task2.bullet_num

    return new_missile


def copy_missile_value(new_missile: Missile, missile: Missile):
    new_missile.location[0:3] = missile.location[0:3]
    new_missile.s_range = missile.s_range
    new_missile.commu_prob = missile.commu_prob
    new_missile.task1.number = missile.task1.number
    new_missile.task1.location[0:3] = missile.task1.location[0:3]
    new_missile.task1.value = missile.task1.value
    new_missile.task1.bullet_num = missile.task1.bullet_num
    new_missile.task2.number = missile.task2.number
    new_missile.task2.location[0:3] = missile.task2.location[0:3]
    new_missile.task2.value = missile.task2.value
    new_missile.task2.bullet_num = missile.task2.bullet_num

#输入两个目标，把B的任务给A
def copy_target(new_target: Target, target: Target):
    new_target.number = target.number
    new_target.location[0] = target.location[0]
    new_target.location[1] = target.location[1]
    new_target.value = target.value  # 价值
    new_target.health = target.health  # 当前状态 0 损毁， 1 健康
    new_target.bullet_num = target.bullet_num  # 最小饱和打击量
    new_target.bullet_real = target.bullet_real  # 真实剩余最小饱和打击量


# def dis(a, b):
#     sum2 = 0
#     for i in range(len(a)):
#         sum2 += (a[i]-b[i])**2
#     return sqrt(sum2)


def dis_lon_lat(lon1, lat1, lon2, lat2):
    # lon经度
    # lat纬度
    R = 6371  # 地球半径（千米）
    # 度转弧度
    lon1rad, lat1rad, lon2rad, lat2rad = lon1*pi/180, lat1*pi/180, lon2*pi/180, lat2*pi/180
    d = R * acos(cos(lat1rad) * cos(lat2rad) * cos(lon1rad - lon2rad) + sin(lat1rad) * sin(lat2rad))
    # 百公里误差数米，千公里误差数十米
    return d


def rad_normol(theta):
    if theta > pi:
        theta -= 2*pi
    elif theta < -pi:
        theta += 2*pi
    return theta


def advantage(dan: Missile, ship: Target, hit_theta=None):
    print('修改后的优势值')
    # 距离优势
    if hit_theta is None:  # 如果没有指定打击角度，则生成打击角度
        hit_theta = ship.theta_domain0[0] + (ship.theta_domain0[1]-ship.theta_domain0[0])*random.random()
        hit_theta_abs = rad_normol(ship.theta+hit_theta)  # 绝对角度
        distance = dubins_min_len(dan.location, dan.theta, ship.location, hit_theta_abs, dan.radius)

        hit_theta_tmp = ship.theta_domain1[0] + (ship.theta_domain1[1] - ship.theta_domain1[0]) * random.random()
        hit_theta_abs_tmp = rad_normol(ship.theta+hit_theta_tmp)
        distance_tmp = dubins_min_len(dan.location, dan.theta, ship.location, hit_theta_abs_tmp, dan.radius)
        if distance_tmp < distance:
            hit_theta = hit_theta_tmp
            distance = distance_tmp
        hit_theta_tmp = ship.theta_domain2[0] + (ship.theta_domain2[1] - ship.theta_domain2[0]) * random.random()
        hit_theta_abs_tmp = rad_normol(ship.theta+hit_theta_tmp)
        distance_tmp = dubins_min_len(dan.location, dan.theta, ship.location, hit_theta_abs_tmp, dan.radius)
        if distance_tmp < distance:
            hit_theta = hit_theta_tmp
            distance = distance_tmp
        hit_theta_tmp = ship.theta_domain3[0] + (ship.theta_domain3[1] - ship.theta_domain3[0]) * random.random()
        hit_theta_abs_tmp = rad_normol(ship.theta+hit_theta_tmp)
        distance_tmp = dubins_min_len(dan.location, dan.theta, ship.location, hit_theta_abs_tmp, dan.radius)
        if distance_tmp < distance:
            hit_theta = hit_theta_tmp
            distance = distance_tmp
    else:
        hit_theta_abs = rad_normol(ship.theta + hit_theta)  # 绝对角度
        distance = dubins_min_len(dan.location, dan.theta, ship.location, hit_theta_abs, dan.radius)
        # d_bar = dan.location - ship.location
    # distance = sqrt(d_bar[0] ** 2 + d_bar[1] ** 2)
    if distance > dan.s_range:
        return 0, 0, 0
    else:
        vd = 1 + 500 / distance
    # elif dan.s_range < 50:
    #     vd = 1 + 500 / distance
    # else:
    #     vd = 1 + distance / dan.s_range + dan.s_range / distance

    # 代价优势
    vv = ship.value / dan.value
    value = vd * vv

    return value, hit_theta, distance


def dan_update_task(dan: Missile, target_map):
    # 根据目标实际情况更新导弹上的目标信息

    for i in list(target_map.keys()):
        #加上的
        #为什么要在30这个位置加，并且为什么30再判断任务是否规0
        if dis(dan.location[0:2], target_map[i].location[0:2]) <= 200:  # 更新真实目标情况
            if dan.task1.number:
                if i == dan.task1.number:
                    if target_map[i].health:
                        dan.task1.location[0:2] = target_map[i].location[0:2]
                        dan.task1.bullet_real = target_map[i].bullet_real
                    else:  # 目标已被击毁，取消任务，等待其他任务
                        dan.task1.number = 0

            else:
                if target_map[i].health and (target_map[i].number > 8 or dan.s_range < 100):
                    # 当前无执行目标，又侦察到新目标或者燃料不足
                    val, hit_theta, hit_distance = advantage(dan, target_map[i])
                    if val > 0:
                        # 当前无执行目标，又侦察到新目标或者燃料不足，且满足打击条件，加入打击任务
                        copy_target(dan.task1, target_map[i])
                        dan.advantage_val = val
                        dan.hit_theta = hit_theta
                        dan.hit_distance = hit_distance
                        dan.task_reset()
            if i == dan.task2.number:
                if target_map[i].health:
                    dan.task2.location[0:2] = target_map[i].location[0:2]
                    dan.task2.bullet_real = target_map[i].bullet_real
                else:  # 目标已被击毁，取消任务，等待其他任务
                    dan.task2.number = 0
    num = dan.task1.number
    if num:
        if dis(dan.location[0:2], dan.task1.location[0:2]) < 15 and \
                dis(dan.location[0:2], target_map[num].location[0:2]) > 30:
            # dan丢失了目标，取消任务，等待其他任务
            dan.task1.number = 0

def dan_update_hit(dan: Missile, target_map):
    # 根据目标实际情况更新导弹上的目标信息
    # for i in list(target_map.keys()):
    #     if dis(dan.location[0:2], target_map[i].location[0:2]) <= 2:
    #         # 击中目标
    #         dan.health = 0
    #         dan.hit = 1
    if dan.task1.number and dis(dan.location[0:2], target_map[dan.task1.number].location[0:2]) <= 10:
        # 击中目标
        dan.health = 0
        dan.hit = 1
def intersect(dan1,dan2):
    x1, y1,x4 ,y4 =dan1.location[0],dan1.location[1],dan1.task1.location[0],dan1.task1.location[1]
    x3, y3, x2, y2=dan2.location[0],dan2.location[1],dan2.task1.location[0],dan2.task1.location[1]
    def cross_product(x1, y1, x2, y2):
        return x1 * y2 - y1 * x2

    if (cross_product(x2-x1, y2-y1, x3-x1, y3-y1) * cross_product(x2-x1, y2-y1, x4-x1, y4-y1) < 0) and \
            (cross_product(x4-x3, y4-y3, x1-x3, y1-y3) * cross_product(x4-x3, y4-y3, x2-x3, y2-y3) < 0):
        return False
    else:
        return True
def exchange_decision2(dan1: Missile, dan2: Missile, time):
    # if dan1.number==1 and dan2.number==2:
    #     print('dan1',dan1.task1.number,dan1.advantage_val)
    #     print('dan2', dan2.task1.number, dan2.advantage_val)
    #     print('dan1打2',advantage(dan1, dan2.task1)[0])
    #     print('dan2打1', advantage(dan2, dan1.task1)[0])
    #     print('交换需求判断',dan1.exchange_stage )
    # 交换需求判断

    if  dan1.task1.number != dan2.task1.number :  # 有交换的需求
        copy_target(dan1.exchange_target, dan2.task1)
        copy_target(dan2.exchange_target, dan1.task1)
        dan1.exchange_number = dan2.number

        # if dan1.task1.number == 0:
        #     val1, hit1, distance1 = advantage(dan1, dan2.task1)
        #     print(distance1)
        #     dan2.advantage_val = advantage(dan2, dan2.task1)[0]
        #     if val1 > dan2.advantage_val:
        #         print('交换了1')
        #         # 交换
        #         #  任务类型，阶段，0 交换 dan2.task1.number，弹1优势值，弹2优势值，弹1打击角，弹2打击角,弹1距离，弹2距离
        #         copy_target(dan1.task1, dan1.exchange_target)
        #         copy_target(dan2.task1, dan2.exchange_target)
        #         dan1.task1.number=dan1.exchange_target.number
        #         dan2.task1.number = 0
        #
        # elif dan2.task1.number == 0:
        #     val2, hit2, distance2 = advantage(dan2, dan1.task1)
        #     dan1.advantage_val = advantage(dan1, dan1.task1)[0]
        #     if dan1.advantage_val < val2:
        #         # 交换
        #         copy_target(dan1.task1, dan1.exchange_target)
        #         copy_target(dan2.task1, dan2.exchange_target)
        #         dan2.task1.number=dan2.exchange_target.number
        #         dan1.task1.number = 0
        XJ = intersect(dan1, dan2)
        if dan1.task1.number !=0 and dan2.task1.number !=0 and XJ :
            val12, hit12, distance12 = advantage(dan1, dan2.task1)
            val21, hit21, distance21 = advantage(dan2, dan1.task1)
            dan1.advantage_val=advantage(dan1, dan1.task1)[0]
            dan2.advantage_val = advantage(dan2, dan2.task1)[0]
            vallist=[dan1.advantage_val,val12,dan2.advantage_val,val21]
            index_val=vallist.index(max(vallist))
            # if (val12+val21)>(dan1.advantage_val+dan2.advantage_val):
            if index_val % 2==1:
            # if dan1.advantage_val<val12 and dan2.advantage_val<val21:
                # if dan1.number == 1 and dan2.number == 2:
                # 交换
                text = [2, 1, dan1.task1.number, dan2.task1.number, val12, val21, hit12, hit21, distance12, distance21]
                dan1.exchange_text = text

                copy_target(dan1.task1, dan1.exchange_target)
                copy_target(dan2.task1, dan2.exchange_target)
                dan1.advantage_val = val12
                dan1.hit_theta = hit12
                dan1.hit_distance = distance12

                dan2.advantage_val = val21
                dan2.hit_theta = hit21
                dan2.hit_distance = distance21

                #  任务类型，阶段，dan1.task1.number 交换 dan2.task1.number，弹1优势值，弹2优势值，弹1打击角，弹2打击角,弹1距离，弹2距离
                return Message(dan2.number, dan1.number, time, text)  # 生成消息
    return 0
#  # 1：请求交换，2：回复请求，3：确认执行
# 如果有交换需求，并且两弹任务一不一样，，把dan2的任务1给dan1
def exchange_decision(dan1: Missile, dan2: Missile, time):

    if dan1.exchange_stage == 0 and dan1.task1.number != dan2.task1.number:  # 有交换的需求
        copy_target(dan1.exchange_target, dan2.task1)
        copy_target(dan2.exchange_target, dan1.task1)
        dan1.exchange_flag = 1
        dan1.exchange_number = dan2.number

        if dan1.task1.number == 0:
            val1, hit1, distance1 = advantage(dan1, dan2.task1)
            if val1 > dan2.advantage_val:
                # 交换
                text = [2, 1, 0, dan2.task1.number, val1, 0, hit1, 0, distance1, 0]
                dan1.exchange_text = text
                #  任务类型，阶段，0 交换 dan2.task1.number，弹1优势值，弹2优势值，弹1打击角，弹2打击角,弹1距离，弹2距离
                return Message(dan2.number, dan1.number, time, text)  # 生成消息
        elif dan2.task1.number == 0:
            val2, hit2, distance2 = advantage(dan2, dan1.task1)
            if dan1.advantage_val < val2:
                # 交换
                text = [2, 1, dan1.task1.number, 0, 0, val2, 0, hit2, 0, distance2]
                dan1.exchange_text = text
                #  任务类型，阶段，dan1.task1.number 交换 0，弹1优势值，弹2优势值，弹1打击角，弹2打击角,弹1距离，弹2距离
                return Message(dan2.number, dan1.number, time, text)  # 生成消息
        else:
            val12, hit12, distance12 = advantage(dan1, dan2.task1)
            val21, hit21, distance21 = advantage(dan2, dan1.task1)
            if dan1.advantage_val < val12 and val21 > dan2.advantage_val:
                # 交换
                text = [2, 1, dan1.task1.number, dan2.task1.number, val12, val21, hit12, hit21, distance12, distance21]
                dan1.exchange_text = text
                #  任务类型，阶段，dan1.task1.number 交换 dan2.task1.number，弹1优势值，弹2优势值，弹1打击角，弹2打击角,弹1距离，弹2距离
                return Message(dan2.number, dan1.number, time, text)  # 生成消息
    return 0


def receive_message(dan1: Missile, message: Message, time):
    if message.sender not in dan1.state_time or message.send_time >= dan1.state_time[message.sender]:
        # 判断是否为第一次收到消息，判断收到的信息是否为旧消息，旧消息不会更新连接时间
        dan1.state_time[message.sender] = message.send_time  # 更新连接时间

    #如果发送的是导弹信息，那么更新信息，并且检测是否有共同打击目标，如果合并信息，并更新距离，更新打击角度
    if not isinstance(message.text, list):
        if message.sender not in dan1.connect_time or message.send_time >= dan1.connect_time[message.sender]:
            # 判断是否为第一次收到消息，判断收到的信息是否为旧消息，旧消息会被舍弃,
            dan1.connect_time[message.sender] = message.send_time  # 更新连接时间
            #没有就添加，有就更新
            if message.sender not in dan1.connect_buffer:
                dan2 = copy_missile(message.text)
                dan1.connect_buffer[message.sender] = dan2
            else:
                copy_missile_value(dan1.connect_buffer[message.sender], message.text)
            if message.text.task1.number:
                # 更新相邻打击角信息  # 添加到二维字典中
                add2d_dict(dan1.hit_theta_dict, message.text.task1.number, message.sender, [message.text.hit_theta, time])
            if dan1.task1.number and message.text.task1.number == dan1.task1.number:
                # 更新友弹航程信息
                s = dubins_min_len(message.text.location, message.text.theta, dan1.task1.location, 0, dan1.radius)
                val = message.text.s_range - s
                add2d_dict(dan1.task_value_dict, dan1.task1.number, message.sender, [val, time])  # 添加到二维字典中

                #检测接受和发送导弹是否有同样的任务，如果有，就合并信息，并且更新同伴信息
                if dan1.task1.number in message.text.task_value_dict:
                    for k in message.text.task_value_dict[dan1.task1.number]:
                        if k != dan1.number:
                            if k in dan1.task_value_dict[dan1.task1.number]:
                                # 以剩余航程小的更新
                                if dan1.task_value_dict[dan1.task1.number][k][0] < \
                                        message.text.task_value_dict[dan1.task1.number][k][0]:
                                    dan1.task_value_dict[dan1.task1.number][k] = \
                                        dan1.task_value_dict[dan1.task1.number][k]
                                else:
                                    dan1.task_value_dict[dan1.task1.number][k] = \
                                        message.text.task_value_dict[dan1.task1.number][k]
                            else:
                                dan1.task_value_dict[dan1.task1.number][k] = \
                                    message.text.task_value_dict[dan1.task1.number][k]
                # 更新相邻打击角信息
                if message.sender == dan1.hit_theta_left[1]:
                    dan1.hit_theta_left = [message.text.hit_theta, message.sender, dan1.task1.number]
                elif message.sender == dan1.hit_theta_right[1]:
                    dan1.hit_theta_right = [message.text.hit_theta, message.sender, dan1.task1.number]
                else:
                    if dan1.hit_theta_left[0] < message.text.hit_theta <= dan1.hit_theta:
                        dan1.hit_theta_left = [message.text.hit_theta, message.sender, dan1.task1.number]
                    elif dan1.hit_theta <= message.text.hit_theta < dan1.hit_theta_right[0]:
                        dan1.hit_theta_right = [message.text.hit_theta, message.sender, dan1.task1.number]

            return 0

    #如果是其他列表信息，又可能是广播也有可能是交换，1是买卖也是广播，2是交换，3是联盟
    else:
        if message.text[0] == 1:  # 买卖协议
            if message.text[1] == 1 and dan1.buy_stage == 0 and dan1.buy_need:
                # 第一阶段,广播任务投标，且投标方(buyer)有买入需求，处于空闲,则接收任务到缓存器，随后选择最优的投标
                new_target = Target(0)
                copy_target(new_target, message.text[2])
                if dan1.task2.number == 0:
                    dan1.buy_task_value[message.sender] = (0, 0, 0)
                if dan1.task1.number == 0:
                    dan1.buy_task_value[message.sender] = advantage(dan1, new_target)
                dan1.buy_task_buffer[message.sender] = new_target
                return 0
            elif message.text[1] == 2 and dan1.sell_stage == 0 and dan1.sell_need:  # 第二阶段,评标,且 拍卖方处于第一阶段
                if message.text[2]:  # 执行任务买入，立即回复
                    task1_flag, task2_flag = 0, 0
                    if message.text[2] == dan1.task2.number and dan1.sell_bullet_num > 0:
                        task1_flag = dan1.task2.number
                        dan1.sell_bullet_num -= 1
                    if message.text[3] == dan1.task2.number and dan1.sell_bullet_num > 0:
                        task2_flag = dan1.task2.number
                        dan1.sell_bullet_num -= 1

                    if dan1.sell_bullet_num < 1:  # 拍卖结束进入第三阶段
                        dan1.sell_stage = 3  # 投标进入第3阶段
                        dan1.sell_time = 0
                    if task1_flag or task2_flag:
                        dan1.sell_num += 1
                        dan1.sell_reply_buffer.add(message.sender)  # 等待此对象回复
                        text = [1, 3, task1_flag, task2_flag]
                        return Message(message.sender, dan1.number, time, text)
                    else:
                        return 0
                else:  # 传播任务买入，稍后处理
                    dan1.buy_bidder_buffer[message.sender] = message
                    return 0

            elif message.sender == dan1.buy_number and dan1.buy_stage == 2 and message.text[1] == 3:
                # 第3阶段,执行,且 投标方处于第2阶段
                if message.sender in dan1.buy_reply_buffer:
                    dan1.buy_reply_buffer.discard(dan1.buy_number)  # 从等待中清空
                else:
                    return 0
                if message.text[2] == 0 and message.text[3] == 0:
                    # 未中标
                    dan1.buy_reset()
                    # dan1.buy_stage = 0  # 投标进入第0阶段
                    # dan1.buy_time = 0
                    # dan1.buy_flag = 0
                    # dan1.buy_target.number = 0
                    # dan1.buy_number = 0
                    return 0

                task1_flag, task2_flag = 0, 0
                if message.text[2] == dan1.buy_target.number:
                    task1_flag = message.text[2]
                    copy_target(dan1.task1, dan1.buy_target)
                    dan1.task1.bullet_num = 1
                    dan1.task_reset()
                if message.text[3] == dan1.buy_target.number:
                    task2_flag = message.text[3]
                    copy_target(dan1.task2, dan1.buy_target)
                    dan1.task2.bullet_num = 1
                text = [1, 4, task1_flag, task2_flag]  # 1, 阶段4：中标结果确认， task1成交任务编号，task2成交任务编号

                if task1_flag or task2_flag:
                    dan1.buy_reset()
                    # dan1.buy_stage = 0  # 投标进入第0阶段
                    # dan1.buy_time = 0
                    # dan1.buy_flag = 0
                    # dan1.buy_target.number = 0
                    # dan1.buy_number = 0
                    # 重发标识开启
                    dan1.buy_reply_flag = 1
                    dan1.buy_reply_message = Message(message.sender, dan1.number, time, text)
                    return dan1.buy_reply_message
                else:  # 信息不匹配，继续等待
                    return 0
            elif message.text[1] == 4:
                if message.sender in dan1.sell_reply_buffer:
                    dan1.sell_reply_buffer.discard(message.sender)
                    dan1.sell_num -= 1
                else:
                    return 0
                if message.text[2] == dan1.task2.number:
                    dan1.task2.bullet_num -= 1
                if message.text[3] == dan1.task2.number:
                    dan1.task2.bullet_num -= 1
                if dan1.task2.bullet_num < 1:  # 任务出售完毕
                    dan1.task2.number = 0
                    dan1.sell_reset()
                    # dan1.sell_flag = 0
                    # dan1.sell_stage = 0
                    # dan1.sell_time = 0
                if dan1.sell_num == 0:  # 等待任务确认完毕
                    dan1.sell_reset()
                    # dan1.sell_flag = 0
                    # dan1.sell_stage = 0
                    # dan1.sell_time = 0
                return 0

        elif message.text[0] == 2:  # 交换协议
            if message.text[1] == 1 and dan1.exchange_stage == 0:
                if message.text[3] == dan1.task1.number:
                    dan1.exchange_text = message.text
                    dan1.exchange_text[1] = 2
                    dan1.exchange_stage = 2
                    dan1.exchange_time = 0
                    dan1.exchange_flag = 1
                    dan1.exchange_number = message.sender
                    dan1.exchange_reply_buffer.add(message.sender)
                    text = [2, 2, message.text[2], message.text[3]]
                    return Message(message.sender, dan1.number, time, text)
                else:
                    text = [2, 2, 0, 0]
                    return Message(message.sender, dan1.number, time, text)
            elif message.sender == dan1.exchange_number and message.text[1] == 2 and dan1.exchange_stage == 1:
                dan1.exchange_reset()
                # dan1.exchange_number = 0
                # dan1.exchange_stage = 0
                # dan1.exchange_time = 0
                # dan1.exchange_flag = 0
                if message.text[3] == dan1.exchange_target.number:
                    copy_target(dan1.task1, dan1.exchange_target)
                    dan1.advantage_val = dan1.exchange_text[4]
                    dan1.hit_theta = dan1.exchange_text[6]
                    dan1.hit_distance = dan1.exchange_text[8]
                    dan1.task_reset()
                    text = message.text
                    text[1] = 3
                    dan1.exchange_target.number = 0
                    # 开启重发
                    dan1.exchange_reply_flag = 1
                    dan1.exchange_reply_message = Message(message.sender, dan1.number, time, text)
                    return dan1.exchange_reply_message
                else:
                    return 0
            elif message.sender == dan1.exchange_number and message.text[1] == 3 and dan1.exchange_stage == 2:
                if message.sender in dan1.exchange_reply_buffer:
                    dan1.exchange_reply_buffer.discardt(message.sender)
                else:
                    return 0
                if message.text[3] == dan1.task1.number:
                    copy_target(dan1.task1, dan1.exchange_target)  # 完成交换
                    dan1.advantage_val = dan1.exchange_text[5]
                    dan1.hit_theta = dan1.exchange_text[7]
                    dan1.hit_distance = dan1.exchange_text[9]
                    dan1.task_reset()
                else:
                    # 交换错误，将对方提供的任务放进目标缓存器中
                    new_target = Target(0)
                    copy_target(new_target, dan1.exchange_target)
                    dan1.target_buffer.append(new_target)
                dan1.exchange_number = 0
                # dan1.exchange_stage = 0
                # dan1.exchange_time = 0
                # dan1.exchange_flag = 0
                # dan1.exchange_number = 0
                dan1.exchange_target.number = 0

                return 0

        elif message.text[0] == 3:  # 联盟协议
            if message.text[1] == 1 and message.text[2] == 1 and dan1.ally_stage == 0:  # 请求阶段
                if dan1.ally_num < 4:
                    # 联盟成员数量少于4
                    # message.text[2]标识符为1
                    dan1.ally_tmp_number = message.sender
                    dan1.ally_stage = 2  # 进入等待确认阶段
                    dan1.ally_flag = 1  # 开启计时
                    dan1.ally_reply_buffer.add(message.sender)
                    text = [3, 2, 1]  # 标识符为1同意加入
                else:
                    text = [3, 2, 0]  # 标识符为0拒绝加入
                return Message(message.sender, dan1.number, time, text)
            elif message.text[1] == 2 and dan1.ally_stage == 1:  # 回复阶段
                if message.text[2] == 1:  # 对方同意加入
                    if message.sender == dan1.ally_tmp_number:  # 加入联盟
                        dan1.ally_buffer[message.sender] = 0
                        dan1.ally_num += 1
                        # 状态归零
                        dan1.ally_reset()
                        # dan1.ally_tmp_number = 0
                        # dan1.ally_flag = 0
                        # dan1.ally_time = 0
                        # dan1.ally_stage = 0
                        text = [3, 3, 1]
                        # 开启重发
                        dan1.ally_reply_flag = 1
                        dan1.ally_reply_message = Message(message.sender, dan1.number, time, text)
                        return Message(message.sender, dan1.number, time, text)
                    else:
                        text = [3, 3, 0]  # 拒绝
                        return Message(message.sender, dan1.number, time, text)
                elif message.text[2] == 0 and message.sender == dan1.ally_tmp_number:  # 对方拒绝
                    # 状态归零
                    dan1.ally_reset()
                    # dan1.ally_tmp_number = 0
                    # dan1.ally_flag = 0
                    # dan1.ally_time = 0
                    # dan1.ally_stage = 0
                    return 0
            elif message.text[1] == 3 and dan1.ally_stage == 2 and message.sender == dan1.ally_tmp_number:  # 确认阶段
                if message.sender in dan1.ally_reply_buffer:
                    dan1.ally_reply_buffer.discard(message.sender)
                else:
                    return 0
                if message.text[2] == 1:  # 加入联盟
                    dan1.ally_buffer[message.sender] = 0
                    dan1.ally_num += 1
                # 状态归零
                dan1.ally_reset()
                # dan1.ally_tmp_number = 0
                # dan1.ally_flag = 0
                # dan1.ally_time = 0
                # dan1.ally_stage = 0
                return 0
            elif message.text[1] == 0 and message.text[2] == 1:
                return Message(message.sender, dan1.number, time, dan1)  # 发送自身状态
            elif message.text[1] == 1 and message.text[2] == 0:  # 脱离联盟
                if message.sender in dan1.ally_buffer:
                    dan1.ally_buffer.pop(message.sender)
                    dan1.ally_num -= 1
                    text = [3, 1, 0]
                    return Message(message.sender, dan1.number, time, text)  # 确认脱离
                return 0


################################cloop_data
def millerToXY (lon, lat):
    xy_coordinate = []
    #地球周长
    L = 6381372 *math.pi*2
    #平面展开，将周长视为X轴
    W = L
    #Y轴约等于周长一般
    H = L/2
    #米勒投影中的一个常数，范围大约在正负2.3之间
    mill = 2.3
    #将经度从度数转换为弧度
    x = lon*math.pi/180
    # 将纬度从度数转换为弧度
    y = lat*math.pi/180
    #这里是米勒投影的转换
    y = 1.25*math.log(math.tan(0.25*math.pi+0.4*y))
    # 这里将弧度转为实际距离 ，转换结果的单位是公里
    x = (W/2)+(W/(2*math.pi))*x/1000
    y = (H/2)-(H/(2*mill))*y/1000
    # xy_coordinate.append((int(round(x)),int(round(y))))
    return x,y

def run(input_list):
    dan_map = {}
    target_map = {}
    dan_num,ship_num=int(input_list[0]), int(input_list[1])
    ship_live={}
    dan_live = set()
    #misslie的信息
    for i in range(dan_num):
        index,x,y,theta=input_list[4*i+2],input_list[4*i+3],input_list[4*i+4],input_list[4*i+5]
        index=int(index)
        dan_map[index] = Missile(index)
        dan_map[index].location[0], dan_map[index].location[1]= x,y
        #经纬度坐标转换
        #dan_map[index].location[0],dan_map[index].location [1]=millerToXY(dan_map[index].location[0],dan_map[index].location [1])
        dan_map[index].theta = 0
        dan_live.add(index)

     # ship的信息
    for j in range(ship_num):
        index, x,y, theta, bullet_num= input_list[5*j+4*dan_num+ 2], input_list[5*j+4*dan_num + 3],input_list[5*j+4*dan_num + 4]\
                                , input_list[5*j+4*dan_num + 5], input_list[5*j+4*dan_num + 6]
        index, bullet_num=int(index),int(bullet_num)
        target_map[index]=Target(index)
        target_map[index].location[0],target_map[index].location [1] = x,y
        target_map[index].value = bullet_num*100
        ship_live[index] = target_map[index].value
        target_map[index].bullet_num = bullet_num
        target_map[index].bullet_real = bullet_num
        #经纬度坐标转换
        #target_map[index].location[0],target_map[index].location [1]=millerToXY(target_map[index].location[0],target_map[index].location [1])
        target_map[index].theta = pi

    tmp = 0
    dan_map_list = list(dan_map.keys())
    target_map_list = list(target_map.keys())
    ship_live = copy.deepcopy(sorted(ship_live.items(), reverse=True,key=lambda x: x[1]))  # 按照价值排序

    #优先分配高等级目标
    for i in range(len(ship_live)):
        index=ship_live[i][0]
        for j in range(0, int(target_map[index].bullet_num)):
            tmp1=dan_map_list[tmp]
            dan_map[tmp1].task2 = copy.deepcopy(target_map[index])
            dan_map[tmp1].task2.bullet_num = 1
            dan_map[tmp1].task1.location[0] = target_map[index].location[0]
            dan_map[tmp1].task1.location[1] = target_map[index].location[1]
            tmp += 1
            if tmp+1>len(dan_map_list):
                break
    return dan_map,target_map,dan_live,ship_live



def listdata_test():
    list_test = []
    missile = read_data('data/MissileParameters8.xls')
    target = read_data('data/Target3.xls')
    n_dan=8
    n_tar=1
    list_test.append(n_dan)
    list_test.append(n_tar)
    for i in range(n_dan):
        list_test.append(missile[i + 1][0])
        list_test.append(missile[i + 1][1])
        list_test.append(missile[i + 1][2])
        list_test.append(0)

    for j in range(n_tar):
        list_test.append(target[j + 1][0])
        list_test.append(target[j + 1][1])
        list_test.append(target[j + 1][2])
        list_test.append(pi)
        list_test.append(target[j + 1][8])
    return list_test


#每一次得到新的位置信息等
def get_data( dan_map, target_map, dan_live, ship_live,list_data):
    dan_num, ship_num = int(list_data[0]), int(list_data[1])

    danlive_index=set() #统计目前存在的导弹
    shiplive_index = set()  # 统计目前存在的目标
    #更新导弹信息
    for i in range(dan_num):
        index,x,y,theta=list_data[4*i+2],list_data[4*i+3],list_data[4*i+4],list_data[4*i+5]
        index=int(index)
        danlive_index.add(index)


        if index in list(dan_map.keys()): #更新信息
            dan_map[index].location[0], dan_map[index].location[1] = x, y
            # 经纬度坐标转换
            #dan_map[index].location[0],dan_map[index].location [1]=millerToXY(dan_map[index].location[0],dan_map[index].location [1])
            dan_map[index].theta = theta
        else:
            dan_map[index] = Missile(index)
            dan_map[index].location[0], dan_map[index].location[1] = x, y
            # 经纬度坐标转换
            # dan_map[index].location[0],dan_map[index].location [1]=millerToXY(dan_map[index].location[0],dan_map[index].location [1])
            dan_map[index].theta = theta
            dan_live.add(index)
    # 删除不存在的导弹
    for dan_index in list(dan_map.keys()):
        if dan_index not in danlive_index:
            dan_map.pop(dan_index)
            dan_live.remove(dan_index)

     # ship的信息
    for j in range(ship_num):
        index, x,y, theta, bullet_num= list_data[5*j+4*dan_num+ 2], list_data[5*j+4*dan_num + 3],list_data[5*j+4*dan_num + 4]\
                                , list_data[5*j+4*dan_num + 5], list_data[5*j+4*dan_num + 6]
        index, bullet_num=int(index),int(bullet_num)
        shiplive_index.add(index)
        if index in list(target_map.keys()):
            target_map[index].location[0], target_map[index].location[1] = x, y
            target_map[index].value = bullet_num * 100

            target_map[index].bullet_num = bullet_num
            target_map[index].bullet_real = bullet_num
            # 经纬度坐标转换
            # target_map[index].location[0],target_map[index].location [1]=millerToXY(target_map[index].location[0],target_map[index].location [1])
            target_map[index].theta = theta

        else:
            target_map[index]=Target(index)
            target_map[index].location[0], target_map[index].location[1] = x, y
            target_map[index].value = bullet_num * 100

            target_map[index].bullet_num = bullet_num
            target_map[index].bullet_real = bullet_num
            # 经纬度坐标转换
            # target_map[index].location[0],target_map[index].location [1]=millerToXY(target_map[index].location[0],target_map[index].location [1])
            target_map[index].theta = theta
            ship_live.append((index, target_map[index].value))

            #分配任务
            for j in range(0, int(target_map[index].bullet_num)):
                for jj in list(dan_map.keys()):
                    if dan_map[jj].task2.number == 0:
                        dan_map[jj].task2 = copy.deepcopy(target_map[index])
                        dan_map[jj].task2.bullet_num = 1
                        break

    for ship_index in list(target_map.keys()):
        if ship_index not in shiplive_index:
            target_map[ship_index].health=0
    for n_ship_live in range(len(ship_live) - 1, -1, -1):
        num = ship_live[n_ship_live][0]
        if target_map[num].health == 0:
            ship_live.pop(n_ship_live)

    return  dan_map,target_map,dan_live,ship_live

def allocation(dan_map,target_map,dan_live,ship_live,message_global,time):
    dan_map_keys = list(dan_map.keys())
    for a in range(len(dan_map_keys) - 1):
        for b in range(a + 1, len(dan_map_keys)):
            i = dan_map_keys[a]
            j = dan_map_keys[b]
            dis1 = dis(dan_map[i].location, dan_map[j].location)
            dis1=0
            if dis1 < dan_map[i].commu_distance * dan_map[i].commu_prob:  # i 的信息可以传输给j
                if time % 2 == 1:  # 隔1个时隙发送一次导弹状态
                    message = Message(dan_map[j].number, dan_map[i].number, time, dan_map[i])
                    message_global.append(message)

                message = dan_map[i].broadcast(dan_map[j].number, time)  # 广播任务
                if message:
                    message_global.append(message)
                    # flag = isinstance(message.text, list)
                if i < j:  # 编号小的发起交换
                    message = exchange_decision2(dan_map[i], dan_map[j], time)
                    if message:
                        message_global.append(message)
            if dis1 < dan_map[j].commu_distance * dan_map[j].commu_prob:  # j 的信息可以传输给i
                if time % 2 == 1:
                    message = Message(dan_map[i].number, dan_map[j].number, time, dan_map[j])
                    message_global.append(message)

                message = dan_map[j].broadcast(dan_map[i].number, time)  # 广播任务
                if message:
                    message_global.append(message)
                if j < i:  # 编号小的发起交换
                    message = exchange_decision(dan_map[j], dan_map[i], time)
                    if message:
                        message_global.append(message)
    # time_end = time_py.clock()  # 记录结束时间
    # time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    # print('通信时间：', time_sum, 'dan:', len(dan_map))

    # 信息传播
    message_global.sort(key=lambda x: x.receive_time, reverse=False)  # 信息按照送达时间排序

    while message_global and message_global[0].receive_time <= time:
        message = message_global.pop(0)
        if random.random() <=1 :  # 一定概率丢失或者中断   con_para['com_pro']
            if message.receiver in dan_map:  # 接收方还存在
                new_message = receive_message(dan_map[message.receiver], message, time)
                if new_message:
                    message_global.append(new_message)  # 生成决定

    for i in list(dan_map.keys()):
        dan_update_task(dan_map[i], target_map)  # 弹根据观测更新目标状态
        dan_map[i].moving2(1)  # 运动
        # tools.task_apportion(dan_map[i])  # 任务更新

        dan_map[i].ally_request(time)  # 加入联盟请求
        dan_map[i].ally(time)  # 联盟信息维护
        dan_map[i].ageing(time)   # 时效性检测

        message = dan_map[i].reply()  # 三次重发
        if message:
            for mess in message:
                message_global.append(mess)

        message = dan_map[i].buy_doing(time)  # 买卖信息操作
        if message:
            message_global.append(message)
        if dan_map[i].sell_flag:
            message = dan_map[i].bidder(time)
            if message:
                for mess in message:
                    message_global.append(mess)

        dan_update_hit(dan_map[i], target_map)  # 弹根据观测更新目标状态
        if dan_map[i].hit:
            tar_num = dan_map[i].task1.number
            if tar_num:  # tar_num != 0
                target_map[tar_num].bullet_num -= 1
                target_map[tar_num].bullet_real -= 1
                if target_map[tar_num].bullet_real < 1:
                    target_map[tar_num].health = 0
                    # target_map[tar_num].value = 0
                    print("目标：%d 被击毁" % tar_num)

        if dan_map[i].task1.number == 0:
            # 向潜在目标飞行
            for n_ship_live in range(len(ship_live)-1, -1, -1):
                num = ship_live[n_ship_live][0]
                if target_map[num].health:
                    if advantage(dan_map[i], target_map[num])[0] > 0:
                        dan_map[i].task1.location[0] = target_map[num].location[0]
                        dan_map[i].task1.location[1] = target_map[num].location[1]
                        break
                # else:
                #     ship_live.pop(n_ship_live)


    return dan_map,target_map,dan_live,ship_live,message_global

def movingdata(dan_map,target_map,time):
    dt=1
    list_data=[]
    dan_num=len(dan_map)
    for i in list(dan_map.keys()):
        dan_map[i].moving1(dt)
        if not dan_map[i].health:
            dan_num=dan_num-1
    target_num=len(target_map.keys())
    for i in list(target_map.keys()):
        target_map[i].moving(dt)
        if not target_map[i].health:
            target_num=target_num-1

    #模拟新任务出现

    list_data.append(dan_num)
    list_data.append(target_num)
    for i in list(dan_map.keys()):
        if dan_map[i].health:
            list_data.append(dan_map[i].number)
            list_data.append(dan_map[i].location[0])
            list_data.append(dan_map[i].location[1])
            list_data.append(dan_map[i].theta)


    for i in list(target_map.keys()):
        if  target_map[i].health:
            list_data.append(target_map[i].number)
            list_data.append(target_map[i].location[0])
            list_data.append(target_map[i].location[1])
            list_data.append(target_map[i].theta)
            list_data.append(target_map[i].bullet_num)

    if time==80:
        target = read_data('data/Target3.xls')
        list_data[1]=list_data[1]+1
        list_data.append(4)
        list_data.append(target[2][1])
        list_data.append(target[2][2])
        list_data.append(pi)
        list_data.append(target[2][8])

    if time==100:
        target = read_data('data/Target3.xls')
        list_data[1]=list_data[1]+1
        list_data.append(5)
        list_data.append(target[3][1])
        list_data.append(target[3][2])
        list_data.append(pi)
        list_data.append(target[3][8])

    return list_data

def result_allocate(dan_map):
    list_out=[]
    dan_map_list=list(dan_map.keys())
    for i in dan_map_list:
        list_out.append(i)
        list_out.append(dan_map[i].task1.number)
        if dan_map[i].hit_theta < 0:
            dan_map[i].hit_theta += 2 * pi
        list_out.append(dan_map[i].hit_theta)

if __name__ == '__main__':
    con_para = {
        'com_pro': 1,
        'shoot_pro': 1,
    }
    # 初始化数据
    time_slot = 20
    time_end = 3001
    dt = 1  # 飞行器采样频率

    message_global = []  # 全局信息流
    list_test=listdata_test()
    dan_map, target_map, dan_live, ship_live = run(list_test)
    message_global = []  # 全局信息流
    for time in range(1, time_end):
        print(time)
        #模拟每次得到数据，导弹的位置和速度和是否存活，和目标位置
        if time != 1:
            #传输的数据
            list_data=movingdata(dan_map,target_map,time)
            #目标区域 目标位置，
            dan_map, target_map, dan_live, ship_live=get_data(dan_map, target_map, dan_live, ship_live,list_data)
        dan_map, target_map, dan_live, ship_live, message_global=allocation(dan_map,target_map,dan_live,ship_live,message_global,time)
        #分配结果
        list_out=result_allocate(dan_map)
        if len(dan_live) < 1 or len(ship_live) < 1:
            time_end = time
            break

    print('剩余飞行器数量：', len(dan_map))
    for i in list(target_map.keys()):
        if target_map[i].health:
            print("目标：%d 存活, 剩余饱和打击量：%d" % (target_map[i].number, target_map[i].bullet_real))