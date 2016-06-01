#include<Windows.h>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<iomanip>
#include<regex>
using namespace std;
ofstream fout;
ifstream fin;
std::regex CPPFILE(R"A(.*(\.cpp|formula)$)A");

static std::wstring CPToUTF16(unsigned code_page, const std::string& input)
{
    auto const size = MultiByteToWideChar(code_page, 0, input.data(), static_cast<int>(input.size()), nullptr, 0);

    std::wstring output;
    output.resize(size);

    if (size == 0 || size != MultiByteToWideChar(code_page, 0, input.data(), static_cast<int>(input.size()), &output[0], static_cast<int>(output.size())))
        output.clear();

    return output;
}

std::wstring UTF8ToUTF16W(const std::string &input)
{
    return CPToUTF16(CP_UTF8, input);
}
//CP_ACP
std::string UTF16ToUTF8(const std::wstring& input,unsigned codepage = CP_UTF8)
{
    auto const size = WideCharToMultiByte(codepage, 0, input.data(), static_cast<int>(input.size()), nullptr, 0, nullptr, nullptr);

    std::string output;
    output.resize(size);

    if (size == 0 || size != WideCharToMultiByte(codepage, 0, input.data(), static_cast<int>(input.size()), &output[0], static_cast<int>(output.size()), nullptr, nullptr))
        output.clear();

    return output;
}

stringstream ss;

void merge(string root,int d=0)
{
    WIN32_FIND_DATAW fd;
    HANDLE hD = FindFirstFileW(UTF8ToUTF16W(root+"*").c_str(),&fd);
    string FileName;
    do{
        FileName = UTF16ToUTF8(fd.cFileName);
        if( FileName[0] == '.' )continue;
        if( fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
        {
            cout<<setw(d*4)<<""<<FileName<<endl;
            ss<<setw(d*2)<<""<<"- "<<FileName<<'\n';
            merge(root+FileName+"\\",d+1);
            continue;
        }
        else if( regex_match(FileName,CPPFILE) )
        {
            cout<<setw(d*4)<<""<<FileName<<endl;
            ss<<setw(d*2)<<""<<"- "<<FileName<<'\n';
            fin.open(root+UTF16ToUTF8(fd.cFileName,CP_ACP));
            string buf((std::istreambuf_iterator<char>(fin)),
                        std::istreambuf_iterator<char>());
            fin.close();
            fout<<"Code:"<<root<<FileName<<endl
                <<"================"<<endl
                <<endl
                <<"```cpp"<<endl
                <<buf<<endl
                <<"```"<<endl
                <<endl;
        }
    }while(FindNextFileW(hD,&fd));
}
int main()
{
    ss<<"# 2016-NTHU-CodeBook\n\n";
    fout.open("CodeBook.md");
    fout<<"Codebook\n=======\n\n";
    
    merge(".\\");
    
    fout.close();
    cout<<endl;
    fout.open("README.md");
    fout<< ss.str();
    fout.close();
}
