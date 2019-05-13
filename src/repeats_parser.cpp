#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <sstream>
#include <tuple>

#include "PAFObject.cpp"

namespace repeats_parser {

	bool parse(std::vector<std::tuple<std::string, int, int>> &repeats, const std::string& path) {
		std::ifstream input(path);
		if (!input) {
			return false;
		}
		std::string line;
		std::string name;
		while (std::getline(input, line)) {
			std::vector<std::string> result;
			std::istringstream iss(line);
			for(std::string s; iss >> s; ) {
				result.push_back(s);
			}
			if (result.size() != 11){
				return false;
			} 
			name.assign(result[0], 1, result[0].size()-1);
			repeats.emplace_back(name, std::atoi(result[9].c_str()), std::atoi(result[10].c_str()));
		}
		return true;
	}

	void remove_covered(std::vector<std::tuple<std::string, int, int>> &repeats, std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
		auto rpt_cmp = [](const std::tuple<std::string, int, int>& a, const std::tuple<std::string, int, int>& b) { 
			if (std::get<2>(a) == std::get<2>(b)) {
				return std::get<1>(a) > std::get<1>(b);
			}
			return std::get<1>(a) < std::get<1>(b);
		};
		std::sort(repeats.begin(), repeats.end(), rpt_cmp);
		std::vector<std::unique_ptr<PAFObject>>::iterator current_paf = paf_objects.begin();
		std::vector<std::tuple<std::string, int, int>>::iterator current_rpt = repeats.begin();
		while(current_paf != paf_objects.end()) {
			while(std::get<2>(*current_rpt) <= (*current_paf)->t_end) {
				if (std::get<1>(*current_rpt) >= (*current_paf)->t_begin) {
					current_rpt = repeats.erase(current_rpt);
				} else {
					current_rpt++;
				}
				if (current_rpt == repeats.end()) return;
			}
			if (std::get<1>(*current_rpt) >= (*current_paf)->t_begin && std::get<1>(*current_rpt) <= (*current_paf)->t_end) {
				int i = 1;
				do {
					if (current_paf + i == paf_objects.end()) return;
					if ((*(current_paf + i))->t_begin <= (*current_paf)->t_end && std::get<2>(*current_rpt) <= (*(current_paf + i))->t_end) {
						current_rpt = repeats.erase(current_rpt);
						break;
					}
					else {
						i++;
					}
				} while ((*(current_paf + i))->t_begin <= (*current_paf)->t_end);
				
			}
			current_paf++;
		}
	}

	std::tuple<int, int, int, int> check_repeats(std::vector<std::tuple<std::string, int, int>> &repeats, const std::string& reference) {
		int i, j;
		std::string first, second;
		for (i = 0; i < repeats.size()-1; i++) {
			for (j = i + 1; j < repeats.size(); j++) {
				if ((std::get<2>(repeats[i]) - std::get<1>(repeats[i])) != (std::get<2>(repeats[j]) - std::get<1>(repeats[j]))) continue;
				first = reference.substr(std::get<1>(repeats[i]), (std::get<2>(repeats[i]) - std::get<1>(repeats[i])));
				second = reference.substr(std::get<1>(repeats[j]), (std::get<2>(repeats[j]) - std::get<1>(repeats[j])));
				if (first.compare(second) == 0) {
					return std::make_tuple(std::get<1>(repeats[i]), std::get<2>(repeats[i]), std::get<1>(repeats[j]), std::get<2>(repeats[j]));
				}
			}
		}
		return std::make_tuple(-1,-1,-1,-1);
	}

} 