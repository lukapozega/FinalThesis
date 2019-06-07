#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <memory>
#include <algorithm>
#include <sstream>
#include <tuple>

#include "PAFObject.cpp"
#include "FASTAQObject.cpp"

namespace repeats_parser {

	bool parse(std::vector<std::tuple<std::string, int, int>> &repeats, const std::string& path) {
		std::ifstream input(path);
		if (!input) {
			return false;
		}
		std::string line;
		std::string name;
		std::string loc;
		int start, end;
		while (std::getline(input, line)) {
			name = line.substr(1, line.find(" ")-1);
			loc = line.substr(line.find(":") + 1);
			start = std::atoi(loc.substr(0, loc.find("-")).c_str());
			end = std::atoi(loc.substr(loc.find("-") + 1).c_str());
			repeats.emplace_back(name, start, end);
		}
		return true;
	}

	void remove_covered(std::vector<std::tuple<std::string, int, int>> &repeats, std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
		auto rpt_cmp = [](const std::tuple<std::string, int, int>& a, const std::tuple<std::string, int, int>& b) { 
			if (std::get<0>(a) == std::get<0>(b)) {
				if (std::get<1>(a) == std::get<1>(b)) {
					return std::get<2>(a) < std::get<2>(b);
				}
				return std::get<1>(a) < std::get<1>(b);
			}
			return std::get<0>(a) > std::get<0>(b);
		};
		std::sort(repeats.begin(), repeats.end(), rpt_cmp);
		std::vector<std::unique_ptr<PAFObject>>::iterator current_paf = paf_objects.begin();
		std::vector<std::tuple<std::string, int, int>>::iterator current_rpt = repeats.begin();
		std::vector<std::string> ids;
		for (auto const &p : paf_objects) {
			if (std::find(ids.begin(), ids.end(), p->t_name) == ids.end()){
				ids.emplace_back(p->t_name);
			}
		}
		int i;
		for (auto const &name : ids) {
			while (std::get<0>(*current_rpt) != name) {
				current_rpt++;
			}
			while(current_paf != paf_objects.end() && (*current_paf)->t_name == name) {
				while(std::get<2>(*current_rpt) <= (*current_paf)->t_end && std::get<0>(*current_rpt) == (*current_paf)->t_name) {
					if (std::get<1>(*current_rpt) >= (*current_paf)->t_begin) {
						current_rpt = repeats.erase(current_rpt);
					} else {
						current_rpt++;
					}
					if (current_rpt == repeats.end()) return;
				}
				i = 1;
				while (current_paf + i != paf_objects.end() && (*(current_paf + i))->t_begin < (*(current_paf))->t_end && std::get<0>(*current_rpt) == name) {
					if (std::get<1>(*current_rpt) < (*(current_paf))->t_end &&
						std::get<2>(*current_rpt) > (*(current_paf + 1))->t_begin &&
						std::get<2>(*current_rpt) < (*(current_paf + 1))->t_end) {
						current_rpt = repeats.erase(current_rpt);
						continue;
					}
					i++;
				}
				current_paf++;
			}
		}
	}

	int alignment(std::string query, std::string target) {
		int n = query.length();
		int m = target.length();
		int i, j;
		int **matrix = new int*[n+1];
		for (int i = 0 ; i <= n ; i++) matrix[i] = new int[m+1];

		for (i = 0; i <= n; i++) {
			for (j = 0; j <= m; j++) {
				matrix[i][j] = 0;
			}
		}

		for (i = 0; i <= n; i++) {
			matrix[i][0] = -i;
		}
	
		for (j = 1; j <= m; j++) {
			matrix[0][j] = -j;
		}

		for (i = 1; i <= n; i++) { 
	        for (j = 1; j <= m; j++) {
	            if (query[i-1] == target[j-1]) {
	                matrix[i][j] = matrix[i - 1][j - 1]; 
	            } else { 
	                matrix[i][j] = std::max({matrix[i - 1][j - 1] - 1,  
	                                matrix[i - 1][j] - 1,  
	                                matrix[i][j - 1] - 1}); 
	            } 
        	} 
    	}
    	return matrix[n][m];
	}

	void check_repeats(std::vector<std::tuple<std::string, int, int>> &repeats, std::vector<std::unique_ptr<FASTAQEntity>>& ref_objects) {
		int i, j;
		std::string first, second;
		std::vector<std::tuple<std::string, int, int>>::iterator current_rpt = repeats.begin();
		std::vector<std::tuple<std::string, int, int>>::iterator it;
		bool rm = true;
		while (current_rpt != repeats.end()) {
			it = current_rpt;
			rm = true;
			while (it != repeats.end()) {
				if (std::get<0>(*current_rpt) == std::get<0>(*it) && it != current_rpt) {
					for (auto const& ref : ref_objects) {
						if (ref->name == std::get<0>(*current_rpt)) {
							first = ref->sequence.substr(std::get<1>(*current_rpt), (std::get<2>(*current_rpt) - std::get<1>(*current_rpt)));
							second = ref->sequence.substr(std::get<1>(*it), (std::get<2>(*it) - std::get<1>(*it)));
							break;
						}
					}
					if (-alignment(first, second) < 0.15 * std::min(std::get<2>(*current_rpt) - std::get<1>(*current_rpt), (std::get<2>(*it) - std::get<1>(*it)))) {
						rm = false;
					}
				}
				it++;
			}
			if (rm) {
				current_rpt = repeats.erase(current_rpt);
				if (current_rpt == repeats.end()) return;
				continue;
			}
			current_rpt++;
		}
		return;
	}

} 