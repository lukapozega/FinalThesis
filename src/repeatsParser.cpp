#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <sstream>

namespace repeatsParser {

	// template<class T> 
	bool parse(std::vector<std::unique_ptr<T>> &objects, const std::string& path) {
		std::ifstream input(path);
		if (!input) {
			std::cerr << "Cannot open the File : "<< path << std::endl;
			return false;
		}
		std::string line;
		std::name;
		while (std::getline(input, line)) {
			std::vector<std::string> result;
			std::istringstream iss(line);
			for(std::string s; iss >> s; ) {
				std::cout << s << std::endl;
				result.push_back(s);
			}
			name.assign(s, 1, s.size()-1);
			objects.push_back(new T(name, result[9], result[10]))
		}
		return true;
	}

} 