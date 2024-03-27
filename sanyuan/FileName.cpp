#include <iostream>
#include <vector>
#include <cmath> // For math functions
#include <ctime> // For time
#include <cstdlib> // For srand and rand
#include <unordered_map>
#include "libxl.h" // ���� LibXL ͷ�ļ�
using namespace libxl;


// �������
double dis(std::vector<double>& a, std::vector<double>& b) {
    double sum2 = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum2 += pow((a[i] - b[i]), 2);
    }
    return sqrt(sum2);
}

// ��ȡ Excel �ļ������ض�ά�ַ�������
std::vector<std::vector<double>> read_data(const wchar_t* path) {
    Book* book = xlCreateBook(); // ����һ�� Excel �ĵ�����
    std::vector<std::vector<double>> data;

    if (book->load(path)) { // ���� Excel �ļ�
        Sheet* sheet = book->getSheet(0); // ��ȡ��һ��������
        if (sheet) {
            //int rowCount = sheet->lastRow(); // ��ȡ����
            //int colCount = sheet->lastCol(); // ��ȡ����

            int rowCount = sheet->lastRow(); // ��ȡ����
            int colCount = sheet->lastCol(); // ��ȡ����
            for (int i = 0; i <= rowCount; ++i) {
                std::vector<double> rowData;
                for (int j = 0; j <= colCount; ++j) {
                    double value = sheet->readNum(i, j); // ��ȡ��Ԫ������
                    std::cout << "�����"<<i<<"�е�"<<j<<"��" << value << "\t"; // �������
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

    book->release(); // �ͷ� Excel �ĵ�����
    return data;
}

//��ӡ����
/* ʹ�÷�ʽ
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
    // ���� con_para �ֵ�
    std::unordered_map<std::string, double> con_para = {
        {"com_pro", 1},
        {"shoot_pro", 1}
    };

    // ��ʼ������
    int time_slot = 20;
    int time_end = 3001;
    int dt = 1;  // ����������Ƶ��
    std::cout << "con_para[\"com_pro\"] = " << con_para["com_pro"] << std::endl;


    //��ȡ�������ݺ�Ŀ������ dan_list[][3]����ţ�x��y�� missle_list[][4]����ţ�x��y����ֵ��
    //����
    const wchar_t* filename = L"MissileParameters8.xls";
    Book* book = xlCreateBook(); // ����һ�� Excel �ĵ�����
    int n_dan = 8;
    double dan_list[8][3];
    if (book->load(filename)) { // ���� Excel �ļ�
        Sheet* sheet = book->getSheet(0); // ��ȡ��һ��������
        if (sheet) {
            //int rowCount = sheet->lastRow(); // ��ȡ����
            //int colCount = sheet->lastCol(); // ��ȡ����

            for (int i = 0; i < 8; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    dan_list[i][j] = sheet->readNum(i + 1, j); // ��ȡ��Ԫ������
                }

            }
        }
        else {
            std::cerr << "Error accessing sheet!" << std::endl;
        }
    }
     book->release(); // �ͷ� Excel �ĵ�����

     //Ŀ��
     const wchar_t* filename2 = L"Target3.xls";
     Book* book2 = xlCreateBook(); // ����һ�� Excel �ĵ�����
     int n_missle = 3;
     double missle_list[3][4];
     if (book2->load(filename2)) { // ���� Excel �ļ�
         Sheet* sheet = book->getSheet(0); // ��ȡ��һ��������
         if (sheet) {
             //int rowCount = sheet->lastRow(); // ��ȡ����
             //int colCount = sheet->lastCol(); // ��ȡ����

             for (int i = 0; i < 3; ++i)
             {
                 for (int j = 0; j < 3; ++j)
                 {
                     missle_list[i][j] = sheet->readNum(i + 1, j); // ��ȡ��Ԫ������
                 }
                 missle_list[i][3] = sheet->readNum(i + 1, 8);
             }
         }
         else {
             std::cerr << "Error accessing sheet!" << std::endl;
         }
     }
     book2->release(); // �ͷ� Excel �ĵ�����

     int rows = sizeof(missle_list) / sizeof(missle_list[0]);
     int cols = sizeof(missle_list[0]) / sizeof(missle_list[0][0]);

     printArray(missle_list, rows, cols);


    return 0;


}