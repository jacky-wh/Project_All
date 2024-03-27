#include <iostream>
#include <vector>
#include <cmath> // For math functions
#include <ctime> // For time
#include <cstdlib> // For srand and rand
#include <unordered_map>
#include "libxl.h" // 包含 LibXL 头文件
using namespace libxl;


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
    std::cout << "con_para[\"com_pro\"] = " << con_para["com_pro"] << std::endl;


    //读取导弹数据和目标数据 dan_list[][3]（编号，x，y） missle_list[][4]（编号，x，y，价值）
    //导弹
    const wchar_t* filename = L"MissileParameters8.xls";
    Book* book = xlCreateBook(); // 创建一个 Excel 文档对象
    int n_dan = 8;
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

     //目标
     const wchar_t* filename2 = L"Target3.xls";
     Book* book2 = xlCreateBook(); // 创建一个 Excel 文档对象
     int n_missle = 3;
     double missle_list[3][4];
     if (book2->load(filename2)) { // 加载 Excel 文件
         Sheet* sheet = book->getSheet(0); // 获取第一个工作表
         if (sheet) {
             //int rowCount = sheet->lastRow(); // 获取行数
             //int colCount = sheet->lastCol(); // 获取列数

             for (int i = 0; i < 3; ++i)
             {
                 for (int j = 0; j < 3; ++j)
                 {
                     missle_list[i][j] = sheet->readNum(i + 1, j); // 读取单元格数据
                 }
                 missle_list[i][3] = sheet->readNum(i + 1, 8);
             }
         }
         else {
             std::cerr << "Error accessing sheet!" << std::endl;
         }
     }
     book2->release(); // 释放 Excel 文档对象

     int rows = sizeof(missle_list) / sizeof(missle_list[0]);
     int cols = sizeof(missle_list[0]) / sizeof(missle_list[0][0]);

     printArray(missle_list, rows, cols);


    return 0;


}