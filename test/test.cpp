#include<experimental/filesystem>
#include<algorithm>
#include<iostream>
#include<iomanip>
#include<string>
#include<vector>

namespace fs = std::experimental::filesystem::v1;
using namespace std::string_literals;

#define TESTDIR "./test"

#if defined CPP11
    #define CPPV "c++11"
#elif defined CPP14
    #define CPPV "c++14"
#else
    #define CPPV "c++98"
#endif

#define GCCCMD "g++ %s -O2 -std=" CPPV " -o opt.o"

std::vector<const char *> ignore_list {"test",".git","test.cpp"};
struct test_status{
    int total;
    int succ;
    int fail;
    int untest;
};

void ForeachFile(std::string path,test_status &status,int deep = 0)
{
    for(fs::directory_entry data : fs::directory_iterator(path))
    {
        auto p = data.path();
        std::string fname = data.path().filename().string();
        
        if( !fname.empty() && fname[0] == '.')
            continue;
        if( std::find(ignore_list.begin(),ignore_list.end(),fname)!=ignore_list.end() )
            continue;
        if( is_directory(data.path()) )
        {
            std::cout<<std::setw(deep*4)<<""<<p<<std::endl;
            ForeachFile(path + "/" + fname,status,deep+1);
            continue;
        }
        if( p.extension().string()!=".cpp" )
            continue;
        std::cout<<std::setw(deep*4)<<""<<p<<std::endl;
        status.total++;
        
        //find TESTDIR/ path / "test" + fname 
        std::string testCpp = TESTDIR + "/"s + path + "/test" + fname;
        if( fs::exists(testCpp) )
        {
            status.fail++;
        }
        else
        {
            status.untest++;
            std::cout<<std::setw(deep*4+1)<<""<<"test not found"<<std::endl;
        }
    }
}

int main()
{
    std::cout<<"GCCCMD:"<<GCCCMD<<std::endl;
    test_status st {0};
    ForeachFile(".",st);
    std::cout<<"============================="<<std::endl
             <<"total file:"<<st.total<<std::endl
             <<"succ :"<<st.succ<<std::endl
             <<"fail :"<<st.fail<<std::endl
             <<"untest:"<<st.untest<<std::endl
             <<"============================="<<std::endl;
    return st.fail!=0;
}