#include <iostream>
#include <fstream>
#include <string> 
#include <cmath>
#include <ctime>

class Gauss;
class Field;
class Generate;
class Interface;
class Log;
class Data3;
class Slice;
class Wave;
class K_means;
class Control;

//using namespace std;

#define pi 3.14
#define e 2.71
#define eps 1e-6

#define size_field 500
#define size_g_list 100
#define size_data_xy 500
#define size_s_list 30
#define size_k_mens_List 30
#define size_k 50

#pragma pack(push, 1)
struct BMPInfoHeader { //structure for image information
    uint32_t bmpSize = sizeof(BMPInfoHeader); //size of the structure
    int32_t bmpWidth;
    int32_t bmpHeight; 
    uint16_t bmpPlanes = 1;
    uint16_t bmpBitCount = 24; //RGB
    uint32_t bmpCompression = 0; //compression type
    uint32_t bmpSizeImage = 0; 
    int32_t bmpXPelsPerMeter = 0; //horizontal expansion
    int32_t bmpYPelsPerMeter = 0; //vertical expansion
    uint32_t bmpClrUsed = 0; //count of colors used
    uint32_t bmpClrImportant = 0; //count of "important" colors used
}; 

struct BMPHeader { //the structure for the BMP header
    uint16_t bmpType = 0x4D42; //'BM' , defines the file type, 0x4D42 in a 16-bit representation
    uint32_t bmpSize; //size of file in bites
    uint16_t bmpReserved1 = 0; //reservation i = 0
    uint16_t bmpReserved2 = 0; //reservation i = 0
    uint32_t bmpOffBits = sizeof(BMPHeader) + sizeof(BMPInfoHeader); //this is the field that indicates where the bit array starts relative to the beginning of the file, that describes the image
}; 
#pragma pack(pop)

class Gauss{
    public:
    double mu_x;   //Координата по x
    double mu_y;   //Координата по y
    double h;   //Высота
    double sgm_x; //Растяжение по x
    double sgm_y; //Растяжение по y

    Gauss()
    {
        mu_x = 0.0;
        mu_y = 0.0;
        h = 0.0;
        sgm_x = 0.0;
        sgm_y = 0.0;
    }
    Gauss( double g_x,double g_y,double g_h,double g_sigma_x,double g_sigma_y)
    {
        mu_x = g_x;
        mu_y = g_y;
        h = g_h;
        sgm_x = g_sigma_x;
        sgm_y = g_sigma_y;
    }

    void make_gauss( double g_x,double g_y,double g_h,double g_sigma_x,double g_sigma_y)
    {
        mu_x = g_x;
        mu_y = g_y;
        h = g_h;
        sgm_x = g_sigma_x;
        sgm_y = g_sigma_y;
    }

    double value(int i, int j)//Возвращает значениеданной данной Гауссинаны в точке x y
    {
        double intermediate;
        double result;
        double x;
        double y;
        x=i;
        y=j;
        intermediate = -0.5 * ((x - mu_x)*(x - mu_x) / (sgm_x*sgm_x) + (y - mu_y)*(y - mu_y) / (sgm_y*sgm_y));
        if(intermediate < -100 + eps){return eps/2;}
        result = pow(e,intermediate) * (h);
        return result;
    }
    
};
class Field{
    public:
    int size_x;
    int size_y;
    double** Mat = new double*[size_field];

    Field()
    {
        size_x = 0;
        size_y = 0;

        for(int i = 0;i<size_field;++i)
        {
            Mat[i] = new double[size_field];
        }

	for (int i = 0; i < size_field; ++i) {
		for (int j = 0; j < size_field; ++j) {
                        Mat[i][j] = 127.0;
		}
        }
    }

    Field(int h)
    {
        for(int i = 0;i<size_field;++i)
        {
            Mat[i] = new double[size_field];
        }

        for(int i = 0; i < size_field;i++)
        {
            for(int j = 0;j < size_field;j++)
            {
                Mat[i][j] = h;
            }
        }
    }

    Field(int size_xx,int size_yy)
    {
        for(int i = 0;i<size_field;++i)
        {
            Mat[i] = new double[size_field];
        }

        for(int i = 0; i < size_xx;i++)
        {
            for(int j = 0;j < size_yy;j++)
            {
                Mat[i][j] = 127.0;
            }
        }
        size_x = size_xx;
        size_y = size_yy;
    }

    ~Field()
    {
        for(int i = 0;i<size_field;++i)
        {
            delete[] Mat[i];
        }
        delete[] Mat;

    }
    
    void make_field(int size_xx,int size_yy)
    {
        size_x = size_xx;
        size_y = size_yy;
        for(int i = 0; i < size_xx;i++)
        {
            for(int j = 0;j < size_yy;j++)
            {
                Mat[i][j] = 127.0;
            }
        }
    }

    void Fill(Gauss& g)//Добавляет гауссиану в поле
    {
        for(int i = 0;i< size_field;i++)
        {
            for (int j = 0; j < size_field; j++)
            {
                Mat[i][j] = Mat[i][j] + g.value(i,j);
            }
        }
    }

    
};
class Generate{
    public:
    //Gauss G_List[size_g_list];
    Gauss* G_List = new Gauss[size_g_list];
    int len_List;

    Generate()
    {
        Gauss g;
        for (int i = 0; i < size_g_list; i++)
        {
            G_List[i] = g;
        }
        len_List = 0;
    }

    ~Generate()
    {
        delete[] G_List;
    }

    void Add_Gaussian_To_List(const Gauss& g)
    {
        G_List[len_List] = g;
        len_List++;
    }

};
class Config{
    public:
    int Log_interface;
    std::string interface_file_name;
    int Log_server;
    std::string server_file_name;
    std::string GNU_file_name;//
    std::string BMP_file_name;
    //Для поля
    int field_leght;
    int field_width;
    //Для гаусса
    double mu_x;   //Координата по x
    double mu_y;   //Координата по y
    double h;   //Высота
    double sgm_x; //Растяжение по x
    double sgm_y; //Растяжение по y
    //Для Гнуплота
   

    Config()
    {
        std::ifstream f;
        std::string value;
        int i;
        f.open("CONFIG.txt");

        for(i = 1;i < 14;i++)
        {
            f>>value;
            f>>value;
            f>>value;
            

            if(i == 1){
                if(value == "ON"){Log_interface = 1;}
                else{Log_interface = 0;}
            }
            else if(i == 2){interface_file_name = value;}
            else if(i == 3){
                if(value == "ON"){Log_server = 1;}
                else{Log_server = 0;}
            }
            else if(i == 4){server_file_name = value;}
            else if(i == 5){GNU_file_name = value;}
            else if(i == 6){BMP_file_name = value;}
            else if(i == 7){field_leght = std::stoi(value);}
            else if(i == 8){field_width = std::stoi(value);}
            else if(i == 9){mu_x = std::stod(value);}
            else if(i == 10){mu_y = std::stod(value);}
            else if(i == 11){h = std::stod(value);}
            else if(i == 12){sgm_x = std::stod(value);}
            else if(i == 13){sgm_y = std::stod(value);}
            
        }
        f.close();
    }
};
class Log{
    public:
    int Log_Server;// 1 = ON, 0 = OFF
    int Log_Interface;
    std::string server_file_name;
    std::string interface_file_name;

    Log()
    {
        Log_Interface = 0;
        Log_Server = 0;
        server_file_name = "";
        interface_file_name = "";
    }

    Log(int L_S, int L_I, std::string& s_f_n, std::string& i_f_n)
    {
        Log_Interface = L_I;
        Log_Server = L_S;
        server_file_name = s_f_n;
        interface_file_name = i_f_n;
    }

    void log_Time(std::ofstream& out)
    {
        // Получаем текущее время
        std::time_t now = std::time(nullptr);
        // Преобразуем в строку формата "YYYY-MM-DD HH:MM:SS"
        char buffer[100];
        std::strftime(buffer, sizeof(buffer), /*"%Y-%m-%d %H:%M:%S"*/ "%H:%M:%S", std::localtime(&now));    
        // Выводим текущее время
        out << "  " << buffer << std::endl;
    }

    void server(const std::string& message)//Запись в сервер
    {
        if(!Log_Server){return ;}
        std::ofstream out;
        out.open(server_file_name, std::ios::app);
        out << message;
        log_Time(out);
        out.close();
    }

    void interface(const std::string& message)
    {
        if(!Log_Interface){return ;}
        std::ofstream out;
        out.open(interface_file_name, std::ios::app);
        out << message;
        log_Time(out);
        out.close();
    }

    void message_log(const std::string& message)
    {
        server(message);
        interface(message);
    }

    void clean_all()
    {
        std::ofstream out1;//server
        std::ofstream out2;//interface
        out1.open(server_file_name);
        out2.open(interface_file_name);
        out1.close();
        out2.close();
    }
};
class Data3{//Для хранения всех данных для срезов
    public:
    int x_end;
    int y_end;
    //double data_mat[size_data_xy][size_data_xy][3];
    int **A = new int*[size_data_xy];//Подходит ли точка по высоте

    int **B = new int*[size_data_xy];//Не занята ли точка компонент...
    


    Data3()
    {
        x_end = 0;
        y_end = 0;
        for(int i = 0;i < size_data_xy;++i)
        {
            A[i] = new int[size_data_xy];
        }

        for(int i = 0;i < size_data_xy;++i)
        {
            B[i] = new int[size_data_xy];
        }

        for(int i = 0;i<size_data_xy;++i)
        {
            for(int j = 0;j<size_data_xy;++j)
            {
                A[i][j] = 0;
                B[i][j] = 0;
            }
        }
    }


    Data3(Field& field,int height)
    {
        for(int i = 0;i < size_data_xy;++i)
        {
            A[i] = new int[size_data_xy];
        }

        for(int i = 0;i < size_data_xy;++i)
        {
            B[i] = new int[size_data_xy];
        }


        x_end = field.size_x;
        y_end = field.size_y;
        for(int i = 0;i<x_end;++i)
        {
            for(int j = 0;j<y_end;++j)
            {   
                if(field.Mat[i][j] > height+eps){A[i][j] = 1;}
                else{A[i][j] = 0;}
                B[i][j] = 0;
            }
        }
    }

    ~Data3()
    {
        for(int i = 0;i<size_data_xy;++i)
        {
            delete[] A[i];
            delete[] B[i];
        }
        delete[] A;
        delete[] B;
    }
};
class Slice{//Сам срез
    public:
    int x_start;
    int y_start;
    int x_end;
    int y_end;
    int x_mass;
    int y_mass;
    double** cluster = new double*[size_data_xy];
    //double** Mat = new double*[size_field];

    Slice()
    {
        x_start = 0;
        y_start = 0;
        x_end = 0;
        y_end = 0;
        x_mass = 0;
        y_mass = 0;
        for(int i = 0;i<size_data_xy;++i)
        {
            cluster[i] = new double[size_data_xy];
        }
        
        for(int i = 0;i<size_data_xy;++i)
        {
            for(int j = 0;j<size_data_xy;++j){cluster[i][j] = eps/4;}
        }
    }

    ~Slice()
    {
        for(int i = 0;i<size_data_xy;++i)
        {
            delete[] cluster[i];
        }
        delete[] cluster;
    }
    void tuning()
    {   
        int sum_x = 0;
        int sum_y = 0;
        int n = 0;
        for(int i = 0;i<size_data_xy;++i)
        {
            for(int j = 0;j<size_data_xy;++j)
            {
                if(cluster[i][j] > eps)
                {
                    x_start = std::min(x_start,i);
                    y_start = std::min(y_start,j);
                    x_end = std::max(x_end,i);
                    y_end = std::max(y_end,j);
                    sum_x += i;
                    sum_y += j;
                    n++;
                }
            }
        }

        x_mass = sum_x / n;
        y_mass = sum_y / n;
        //std::cout<<x_mass<<" "<<y_mass<<std::endl;
    }
};
class Wave{
    public:
    //Slice List[size_s_list];
    Slice* List = new Slice[size_s_list];
    int list_index;//Длинна списка

    Wave()
    {
        list_index = 0;
    }

    ~Wave()
    {
        delete[] List;
    }
    void wave(Data3& Data,Field& field, int i,int j)//Определяет принадлежность к текущей компонент связности
    {
        if(i < 0 || i >= Data.x_end || j < 0 || j >= Data.y_end){return ;}
        if(Data.A[i][j] == 1  && Data.B[i][j] == 0)
        {
            List[list_index].cluster[i][j] = field.Mat[i][j];
            Data.B[i][j] = 1;
            wave(Data,field,i-1, j-1);
            wave(Data,field,i-1,j);
            wave(Data,field,i - 1,j+1);
            wave(Data,field,i,j+1);
            wave(Data,field,i+1,j+1);
            wave(Data,field,i+1,j);
            wave(Data,field,i+1,j-1);
            wave(Data,field,i,j-1);
        }
        return ;
    }

    void wave_run(Data3& Data,Field& field)
    {   
        list_index = 0;
        for(int i = 0;i<Data.x_end;++i)
        {
            for(int j = 0;j < Data.y_end;++j)
            {
                if(Data.A[i][j] == 1  && Data.B[i][j] == 0)
                {
                    wave(Data,field,i,j);
                    list_index++;
                }
            }
        }
    }

    void tuning()
    {
        for(int i = 0;i<list_index;++i){List[i].tuning();}
    }
    
    void run(Field& field,int height)
    {
        Data3 Data(field,height);
        wave_run(Data,field);
        tuning();
    }
};
class K_means{
    public:
    int array[size_k_mens_List];//Под каждым срезом ставим номер центра, к которому он принадлежит
    int* centers_x = new int[size_k];//Координаты центров
    int* centers_y = new int[size_k];
    int len_k;

    K_means()
    {
        len_k = 0;
    }

    ~K_means()
    {
        delete[] centers_x;
        delete[] centers_y;
    }

    int square_distans(int x1, int y1, int x2, int y2)
    {
        return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y1);
    }

    int find_min(int x,int y,int k)//Находит номер ближайшего центра к данной точке
    {
        int k_min = 0;
        int distance;
        int min_value;
        min_value = square_distans(centers_x[0],centers_y[0],x,y);
        for(int i = 1;i<k;++i)
        {
            distance = square_distans(centers_x[i],centers_y[i],x,y);
            if(distance < min_value)
            {
                k_min = i;
                min_value = distance;
            }
        }
        return k_min;
    }

    void k_means(Slice *List,int len,int k)//Под каждым срезом пишет номер центра,к которому принадлежит
    {
        int n;
        int sum_x;
        int sum_y;
        int flag = 1;

        for(int i = 0;i<k;++i)//Назначили начальные центры
        {
            //centers_x[i] = List[len/k * i + 1].x_mass;
            centers_x[i] = List[i].x_mass;
            centers_y[i] = List[i].y_mass;
        }

        while(flag != k){
            flag = 0;

            for(int i = 0;i<len;++i)//Сопоставили центры
            {
                array[i] = find_min(List[i].x_mass,List[i].y_mass,k);
            }
            /*
            std::cout<<"x"<<array[0]<<" "<<array[1]<<" "<<array[2]<<std::endl;
            std::cout<<"c"<<centers_x[0]<<"/"<<centers_y[0]<<" "<<centers_x[1]<<"/"<<centers_y[1]<<std::endl;
            */
            for(int i = 0;i < k;++i)
            {
                n = 0;
                sum_x = 0;
                sum_y = 0;
                for(int j = 0;j<len;++j)
                {
                    if(array[j] == i)
                    {
                        n++;
                        sum_x += List[j].x_mass;
                        sum_y += List[j].y_mass;
                    }
                }
                if((centers_x[i] == sum_x/n) && (centers_y[i] = sum_y/n)){flag++;}
                else
                {
                    centers_x[i] = sum_x/n;
                    centers_y[i] = sum_y/n;
                }
            }
        }
    }

    void painting(Slice& slice,Field& field,int color)//Раскрашивает поле по срезу в данный цвет(номер)
    {
        for(int i = slice.x_start;i<=slice.x_end;++i)
        {
            for(int j = slice.y_start;j<=slice.y_end;j++)
            {
                if(slice.cluster[i][j] > eps){field.Mat[i][j] = color;}
            }
        }
    }

    void run(Wave& wave,Field& field,int k)
    {
        len_k = k;
        k_means(wave.List,wave.list_index,k);
        for(int i = 0; i < wave.list_index;++i)
        {
            painting(wave.List[i],field,array[i] + 1);
            field.size_x = std::max(field.size_x,wave.List[i].x_end);
            field.size_y = std::max(field.size_y,wave.List[i].y_end);
        }
    }

};
class Control{
	public:
    Log log;
    Config config;
    Generate list;
    Field field;
    Field field_0;//For k_means
    Data3 data;
    Slice slice;
    Wave wave;
    K_means k_m;

	Control(){
        log.Log_Interface = config.Log_interface;
        log.interface_file_name = config.interface_file_name;
        log.Log_Server = config.Log_server;
        log.server_file_name = config.server_file_name;
        field_0 = Field(0);
    }

    void valid_data_Field(int* size_x,int* size_y){
        if(*size_x < eps){
            log.message_log(" Invalid_argument(field)(leght)");
            std::cout<<"Invalid_argument(leght): set by default setting"<<std::endl;
            *size_x = config.field_leght;
        }
        if(*size_y < eps)
        {
            log.message_log(" Invalid_argument(field)(widht)");
            std::cout<<"Invalid_argument(widht): set by default setting"<<std::endl;
            *size_y = config.field_width;
        }
    }

    void valid_data_Gauss(double* sx, double* sy){
    
        //if(*sx < eps || *sy < eps || *h < eps)
        if(*sx < eps)
            {log.message_log(" invalid_argument(gaussian)(sx)");
            std::cout<<"Invalid_argument(sx): set by default\n";
            *sx = config.sgm_x;
        }
        if(*sy < eps)
            {log.message_log(" Invalid_argument(gaussian)(sy)");
            std::cout<<"Invalid_argument(sy): set by default\n";
            *sy = config.sgm_x;
        }
    }

    void INIT(const std::string& method){//INIT
        log.message_log("INIT:");
        if(method == "k")
        {
            int size_x;
            int size_y;

            std::cout << "Enter field dimensions (leght widht): ";
            std::cin >> size_x >> size_y;

            valid_data_Field(&size_x,&size_y);

            field.make_field(size_x,size_y);
            field_0.size_x = size_x;
            field_0.size_y = size_y;
        }

        else if(method == "f")
        {
            std::string file_name;
            int size_x;
            int size_y;
            
            std::cout<<"Enter file name for field data: ";
            std::cin>>file_name;

            std::ifstream file(file_name);

           if(!file){log.message_log(" file opening error(" + file_name + ")");std::cout<<"Can't open "<<file_name<<std::endl;return ;}

            file >> size_x >> size_y;

            valid_data_Field(&size_x,&size_y);

            field.make_field(size_x,size_y);
            field_0.size_x = size_x;
            field_0.size_y = size_y;

            file.close();
        }
        else
        {
            field.make_field(config.field_leght, config.field_width);
            field_0.size_x = config.field_leght;
            field_0.size_y = config.field_width;
        }
        log.message_log(" field created");
    }

    void GS(const std::string& method){//GS
        log.message_log("GS:");
        if(method == "k")
        {
            double x, y, h, sx, sy;
            std::cout << "Enter parameters for Gauss (x y h sx sy ): ";
            
            std::cin >> x >> y >> h >> sx >> sy;

            valid_data_Gauss(&sx,&sy);

            Gauss g(x,y,h,sx,sy);//создаём элемент класса гаусс
            list.Add_Gaussian_To_List(g);//Добавляем его в список
        }
        else if(method == "f")
        {
            std::string file_name;
            double x, y, sx, sy, h;
            Gauss g;
            std::cout<<"Enter name of file for Gauss: ";
            std::cin>>file_name;
            std::ifstream file(file_name);

            if(!file){log.message_log(" file opening error(" + file_name + ")");std::cout<<"Can't open "<<file_name<<std::endl;return ;}

            while(file >> x >> y >> h >> sx >> sy)
            {
                valid_data_Gauss(&sx,&sy);

                g.make_gauss(x,y,h,sx,sy);
                list.Add_Gaussian_To_List(g);
            }
            file.close();
        }
        else
        {
            Gauss g;
            g.make_gauss(config.mu_x, config.mu_y, config.h, config.sgm_x, config.sgm_y);
            list.Add_Gaussian_To_List(g);
        }
        log.message_log(" gaussians added");
    }

    void GNN(){//GNN
        log.message_log("GNN:");
        Gauss g;
        for(int i=0; i < list.len_List;i++)
        {
            g = list.G_List[i];
            field.Fill(g);
        }
        log.message_log(" list of gaussians added into the field");
    }

    void GNU(const std::string& method){//GNU
        log.message_log("GNU:");
        std::string file_name;
        std::ofstream file;
        if(method == "f"){
            std::cout<<"Enter file name for writing field into:";
            std::cin>>file_name;
        }
        else{file_name = config.GNU_file_name;}
        

        file.open(file_name);
        if(!file){log.message_log(" file opening error(" + file_name + ")");std::cout<<"Can't open "<<file_name<<std::endl;return ;}

        for(int i = 0;i < field.size_x;i++)
        {
            for(int j = 0;j < field.size_y;j++)
            {
                file <<i<<" "<<j<<" "<<field.Mat[i][j]<<std::endl;
            }
        }
       log.message_log(" field saved to gnuplot(" + file_name + ")");
        file.close();
    }

    void BMP(const std::string& method){//BMP чёрно-белое
        log.message_log("BMP:");
        std::string file_name;
        int width = field.size_x;
        int height = field.size_y;
        int i;

          if(method == "f"){
            std::cout<<"Enter file name for writing field into:";
            std::cin>>file_name;
        }
        else{file_name = config.BMP_file_name;}

        BMPHeader header;
        BMPInfoHeader infoHeader;

        header.bmpSize = sizeof(BMPHeader) + sizeof(BMPInfoHeader) + width * height * 3;
        infoHeader.bmpWidth = width;
        infoHeader.bmpHeight = -height;

        std::ofstream file(file_name, std::ios::binary);

        file.write(reinterpret_cast<const char*>(&header), sizeof(header));
        file.write(reinterpret_cast<const char*>(&infoHeader), sizeof(infoHeader));

        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
    
                uint8_t red = static_cast<uint8_t>(std::min(std::max(field.Mat[y][x], 0.0), 255.0));// Чтобы не выйти за 256 пикселей
                uint8_t green = static_cast<uint8_t>(std::min(std::max(field.Mat[y][x], 0.0), 255.0));
                uint8_t blue = static_cast<uint8_t>(std::min(std::max(field.Mat[y][x], 0.0), 255.0));

                file.write(reinterpret_cast<const char*>(&blue), sizeof(uint8_t));
                file.write(reinterpret_cast<const char*>(&green), sizeof(uint8_t));
                file.write(reinterpret_cast<const char*>(&red), sizeof(uint8_t));
            }
            for (i = 0; i < (4 - (width * 3) % 4) % 4; i++)
                {
                    uint8_t zero = 0;
                    file.write(reinterpret_cast<const char*>(&zero), sizeof(uint8_t));
                }
        }
        log.message_log(" field saved to bmp(" + file_name + ")");
        file.close();
    }

    void READ_BMP()(Field& field){
        std::string file_name;
        log.message_log("READ BMP:");
        std::cout<<"Enter file name:";
        std::sin>>file_name;
        std::ifstream bmpFile(file_name, std::ios::binary);
        if (!bmpFile) {
            log.message_log(" Error: cannot open file");
            std::cout << "Failed to open BMP file." << std::endl;
            return;
        }

        if (header.bfType != 0x4D42) { // Проверка "BM"
            log.message_log(" Error: file is not bmp");
            std::cout << "File is not bmp!" << std::endl;
            file.close();
            return;
        }
        // Читаем заголовок BMP
        unsigned char header[54];
        bmpFile.read(reinterpret_cast<char*>(header), 54);
        // Получаем ширину и высоту изображения
        int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
        int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);

        // Читаем данные пикселей
        for (int y = height - 1; y >= 0; --y) { // BMP хранит данные снизу вверх
            for (int x = 0; x < width; ++x) {
                unsigned char b = bmpFile.get(); // Читаем B
                unsigned char g = bmpFile.get(); // Читаем G
                unsigned char r = bmpFile.get(); // Читаем R
                double color = (r + g + b) / 3;
                field.Mat[y][x] = value; // Обновляем матрицу значений
            
                bmpFile.ignore((4 - (width * 3) % 4) % 4); // Пропускаем паддинг
            }
            bmpFile.close();
        }
    }

    void DO_SLICE(){
        int height;
        log.message_log("DO_CUT:");
        std::cout<<"Enter the height at which the cut will be made: ";
        std::cin>>height;
        if(height <= 0 || height >= 255){log.message_log(" incorrect data(height)");std::cout<<"height must be between 0 and 255"<<std::endl;}
        log.message_log(" wave start");
        wave.run(field,height);
        log.message_log(" wave end");

    }
    
    void DO_K_MEANS(){
        int k;
        log.message_log("DO_K_MEANS");
        std::cout<<"Enter the number of clusters you want: ";
        std::cin>>k;
        while(k <= 0 || k > wave.list_index)
        {
            log.message_log(" incorrect data(k)");
            std::cout<<"Number must be positive and less then number of slices!"<<std::endl;
            std::cout<<"Try again: ";
            std::cin>>k;
        }
        k_m.run(wave,field_0,k);
    }

    void SAVE_K_MEANS_TO_BMP(){
    
        log.message_log("SAVE_K_MEANS_TO_BMP:");
        std::string file_name;
        int width = field_0.size_y;
        int height = field_0.size_x;
        int k;
        int tmp;

        k = k_m.len_k;
        tmp = (255*255*255)/(k+1);
        std::cout<<"Enter file name for writing field_0 into: ";
        std::cin>>file_name;

        BMPHeader header;
        BMPInfoHeader infoHeader;

        header.bmpSize = sizeof(BMPHeader) + sizeof(BMPInfoHeader) + width * height * 3;
        infoHeader.bmpWidth = width;
        infoHeader.bmpHeight = -height;

        std::ofstream file(file_name, std::ios::binary);

        file.write(reinterpret_cast<const char*>(&header), sizeof(header));
        file.write(reinterpret_cast<const char*>(&infoHeader), sizeof(infoHeader));
        
        for (int x = 0; x < height; ++x)
        {
            for (int y = 0; y < width; ++y)
            {
    
                uint8_t red = static_cast<uint8_t>((tmp * (int)field_0.Mat[x][y]) % 256);// Чтобы не выйти за 256 пикселей
                uint8_t green = static_cast<uint8_t>(((tmp * (int)field_0.Mat[x][y]) / 256) % 256);
                uint8_t blue = static_cast<uint8_t>(((tmp * (int)field_0.Mat[x][y]) / 256) / 256);

                file.write(reinterpret_cast<const char*>(&blue), sizeof(uint8_t));
                file.write(reinterpret_cast<const char*>(&green), sizeof(uint8_t));
                file.write(reinterpret_cast<const char*>(&red), sizeof(uint8_t));
            }
            for (int i = 0; i < (4 - (width * 3) % 4) % 4; i++)
                {
                    uint8_t zero = 0;
                    file.write(reinterpret_cast<const char*>(&zero), sizeof(uint8_t));
                }
        }
        log.message_log(" field_0 saved to bmp(" + file_name + ")");
        file.close();
    }
};
class Interface{
    private:
    Control control;
    public:
    Interface()
    {
        control.log.clean_all();
        control.log.message_log("Interface created");
    }

    void greeting(){
        std::cout << "List of commands:"<<std::endl;
        std::cout << " init( /k/f) - create field \n g(k/f) - add Gaussian \n generate - overlay Gaussians to the field \n gnuplot( /f) - create plot \n bmp( /f) - create BMP file "<<std::endl;
        std::cout
        std::cout << " slice - make a cut of current field" << std::endl;
        std::cout << " k_means - applyind k_means method to current field" << std::endl;
        std::cout << " k_m(bmp) - saving to bmp\n" << std::endl;
        std::cout << " exit - end of the program" << std::endl;
    }

    int execute_command(){
        std::string command;

        std::cout<<"\nEnter command: ";
        std::cin>>command;
        if(command == "init(k)"){control.INIT("k");std::cout<<"Field created"<<std::endl;return 1;}
        else if(command == "init(f)"){control.INIT("f");std::cout<<"Field created"<<std::endl;return 1;}
        else if(command == "init"){control.INIT(" ");std::cout<<"Field created"<<std::endl;return 1;}
        else if(command == "g(k)"){control.GS("k");std::cout<<"Gaussian added to the list"<<std::endl;return 1;}
        else if(command == "g(f)"){control.GS("f");std::cout<<"Gaussian added to the list"<<std::endl;return 1;}
        else if(command == "g"){control.GS(" ");std::cout<<"Gaussian added to the list"<<std::endl;return 1;}
        else if(command == "generate"){control.GNN();std::cout<<"List of gaussians added into the field"<<std::endl;return 1;}
        else if(command == "gnuplot"){control.GNU("");std::cout<<"Gnuplot file created"<<std::endl;return 1;}
        else if(command == "gnuplot()"){control.GNU("");std::cout<<"Gnuplot file created"<<std::endl;return 1;}
        else if(command == "gnuplot(f)"){control.GNU("f");std::cout<<"Gnuplot file created"<<std::endl;return 1;}
        else if(command == "bmp"){control.BMP("");std::cout<<"BMP file created"<<std::endl;return 1;}
        else if(command == "bmp()"){control.BMP("");std::cout<<"BMP file created"<<std::endl;return 1;}
        else if(command == "bmp(f)"){control.BMP("f");std::cout<<"BMP file created"<<std::endl;return 1;}
        else if(command == "exit"){return 0;}
        else if(command == "k_m(bmp)"){control.SAVE_K_MEANS_TO_BMP();std::cout<<"K-Means was saved to bmp"<<std::endl;return 1;}
        else if(command == "k_means"){control.DO_K_MEANS();std::cout<<"Method k_means was applied"<<std::endl;return 1;}
        else if(command == "slice"){control.DO_SLICE();std::cout<<"The cut was made"<<std::endl;return 1;}
        else if(command == "read_bmp"){control.READ_BMP(control.field);return 1;}

        else{
            std::cout<<"This command is not in the list of commands!"<<std::endl;
            execute_command();
        }
        return 1;

    }
    
    void Launch(){

        greeting();

        while(true)
        {
            if(!execute_command()){ control.log.message_log("Exit programm");return ;}
        }
    }
    
};

int main(void)
{
	Interface interface;
    interface.Launch();
    return 0;
}
