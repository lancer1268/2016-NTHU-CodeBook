#include<Windows.h>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<iomanip>
#include<regex>
using namespace std;
ofstream fout;
ifstream fin;
std::regex CPPFILE(R"(.*\.cpp)");

void merge(string root,int d=0)
{
    WIN32_FIND_DATAA fd;
    HANDLE hD = FindFirstFile((root+"*").c_str(),&fd);
    
    do{
        if( fd.cFileName[0] == '.' )continue;
        if( fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
        {
            cout<<setw(d*4)<<""<<fd.cFileName<<endl;
            merge(root+fd.cFileName+"\\",d+1);
            continue;
        }
        else if( regex_match(fd.cFileName,CPPFILE) )
        {
            cout<<setw(d*4)<<""<<fd.cFileName<<endl;
            fin.open(root+fd.cFileName);
            string buf((std::istreambuf_iterator<char>(fin)),
                        std::istreambuf_iterator<char>());
            fin.close();
            fout<<"Code:"<<root<<fd.cFileName<<endl
                <<"================"<<endl
                <<endl
                <<"```cpp"<<endl
                <<buf<<endl
                <<"```"<<endl
                <<endl;
        }
    }while(FindNextFileA(hD,&fd));
}
int main()
{
    fout.open("CodeBook.md");
    fout<<"Codebook\n=======\n\n";
    
    merge(".\\");
    fout.close();
}
