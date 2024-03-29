from math import sin, cos
from xlutils.copy import copy
import xlrd
import xlwt
from numpy import random

# import openpyxl

PI = 3.14152653589793


def produce_data(num):
    # 创建工作workbook
    workbook = xlwt.Workbook()

    # 创建工作表sheet,填入表名Center_error
    worksheet = workbook.add_sheet('UAV_state')

    # 在表1行1列中写入相应的数据11
    for i in range(num):
        worksheet.write(i, 0, 250 + 700 * (random.rand() - 0.5))
        worksheet.write(i, 1, 2800 + 600 * (random.rand() - 0.5))
        worksheet.write(i, 2, 50 + 150 * random.rand())  # 高度
        for j in range(3, 6):  # 速度矢量
            worksheet.write(i, j, 20 * (random.rand() - 0.5))
        for j in range(6, 9):  # 角速度矢量
            worksheet.write(i, j, 0.4 * (random.rand() - 0.5))
        # 四元素
        # 原角度的一半
        theta = (random.rand() - 0.5)
        phi = (random.rand() - 0.5)
        psi = PI * (random.rand() - 0.5)
        q0 = cos(theta) * cos(phi) * cos(psi) + sin(theta) * sin(phi) * sin(psi)
        q1 = sin(phi) * cos(theta) * cos(psi) - cos(phi) * sin(theta) * sin(psi)
        q2 = cos(phi) * sin(theta) * cos(psi) + sin(phi) * cos(theta) * sin(psi)
        q3 = cos(theta) * cos(phi) * sin(psi) - sin(theta) * sin(phi) * cos(psi)
        worksheet.write(i, 9, q0)
        worksheet.write(i, 10, q1)
        worksheet.write(i, 11, q2)
        worksheet.write(i, 12, q3)
        worksheet.write(i, 13, 0)  # 初始代价
        worksheet.write(i, 14, 250 + 500 * (random.rand() - 0.5))
        worksheet.write(i, 15, 2800 + 600 * (random.rand() - 0.5))
        worksheet.write(i, 16, 10 * random.rand())  # 高度

        # if random.rand() > 0.5:
        #     worksheet.write(i, 13, 1)  # 剩余通信量x/1000, 大于1000的为1
        # else:
        #     worksheet.write(i, 13, random.rand())  # 剩余通信量x/1000, 大于1000的为1
        # worksheet.write(i, 14, random.rand())  # 累计的探测成功概率
        # # 与目标的相对位置
        # for j in range(15, 17):
        #     worksheet.write(i, j, 200 * (random.rand() - 0.5))
        # worksheet.write(i, 17, 60 + 100 * random.rand())  # 相对高度

    # 保存表
    workbook.save('./TrainData/Original_data.xlsx')


def read_data(path):
    workbook = xlrd.open_workbook(path)
    table = workbook.sheets()[0]  # 获取第一个sheet表
    row = table.nrows  # 行数
    col = table.ncols  # 列数
    data = []
    for x in range(row):
        data.append(table.row_values(x))
    return data


def produce_controls(path, u, state_t, col):
    try:
        work_book = xlrd.open_workbook(path)
        # table = work_book.sheets()[0]  # 获取第一个sheet表
        # col = table.ncols  # 行数
        new_book = copy(work_book)
        new_table = new_book.get_sheet(0)
        for i in range(4):
            new_table.write(col, i, float(u[i]))
        for i in range(len(state_t)):
            new_table.write(col, i+4, float(state_t[i]))

    except Exception as e:
        # 创建工作workbook
        new_book = xlwt.Workbook()

        # 创建工作表sheet,填入表名Center_error
        new_table = new_book.add_sheet('Controls')
        for i in range(4):
            new_table.write(0, i, float(u[i]))
        for i in range(len(state_t)):
            new_table.write(0, i + 4, float(state_t[i]))
    new_book.save(path)


# def remove_sheet(path):
#     work_book = openpyxl.load_workbook(path)
#     sheet = work_book.active
#     work_book.remove(sheet)
#     work_book.save(path)
#     work_book.close()



if __name__ == '__main__':
    paths = '.\\TrainData\\Controls_data.xlsx'
    #paths = '.\\TrainData\\new.xlsx'
    # u1 = [15, 0.1, 0.1, -0.15]
    # produce_controls(paths, u1)
    produce_data(500)
    # x = read_data('.\\TrainData\\Original_data.xls')
    # remove_sheet(paths)
    # print(x)


