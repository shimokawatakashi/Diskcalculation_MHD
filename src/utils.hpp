#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace Utils 
{

std::vector<std::string> Split(const std::string &s, char delim);

void FileCopy(const std::string filename_in, const std::string filename_out);

bool ConvertStrToBool(std::string strbool);

// template <class T>
class Dictionary
{
public:
    // explicit operator bool() const { return valid; }
    bool Contains(const std::string &key) const { return dict.find(key) != dict.end(); };
    void Insert(std::string key, double item) { dict[key] = item; };

    size_t Size() { return dict.size(); }

    double Get_double(std::string key) 
    { 
        if (Contains(key)) {
            return  dict[key];
        } else {
            return 0.00;
        }
    }

private:
    std::map<std::string, double> dict;
    // bool valid;
};

class InputConfigure
{
public:

    InputConfigure();
    InputConfigure(std::string file_name);
    ~InputConfigure();

    bool Contains(const std::string &key);
    void Insert(std::string key, std::string item);
    void ReadFile(std::string file_name);
    size_t Size();

    std::string Get(std::string key);

    std::string GetInputFileName();

private:

    std::map<std::string, std::string> dict;
    std::string input_filename;

};

} // end namespace Utils

#endif /* UTILS_HPP_ */