#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "Config.h"
#include "bioparser/bioparser.hpp"

class PAFobject {

public:

	std::string query_name;
	std::uint32_t q_length;
	std::uint32_t q_begin;
	std::uint32_t q_end;
	char orientation;
	std::string target_name;
	std::uint32_t t_length;
	std::uint32_t t_begin;
	std::uint32_t t_end;

	PAFobject(
		const char* q_name, std::uint32_t q_name_length,
	    std::uint32_t q_length,
	    std::uint32_t q_begin,
	    std::uint32_t q_end,
	    char orientation,
	    const char* t_name, std::uint32_t t_name_length,
	    std::uint32_t t_length,
	    std::uint32_t t_begin,
	    std::uint32_t t_end,
	    std::uint32_t matching_bases,
	    std::uint32_t overlap_length,
	    std::uint32_t mapping_quality) {
			(this->query_name).assign(q_name, q_name_length);
			this->q_length = q_length;
			this->q_begin = q_begin;
			this->q_end = q_end;
			this->orientation = orientation;
			(this->target_name).assign(t_name, t_name_length);
			this->t_length = t_length;
			this->t_begin = t_begin;
			this->t_end = t_end;
    }

    bool operator < (const PAFobject& object) const {
        if (t_begin == object.t_begin) {
        	return (t_end > object.t_end);
        }
        return (t_begin < object.t_begin);
    }

};

auto paf_cmp = [](const std::unique_ptr<PAFobject>& a, const std::unique_ptr<PAFobject>& b) { 
 		if (a->t_begin == b->t_begin) {
        	return (a->t_end > b->t_end);
        }
        return (a->t_begin < b->t_begin);
	};

void sweepLineAlgorithm(std::vector<std::unique_ptr<PAFobject>> &paf_objects) {
	std::sort(paf_objects.begin(), paf_objects.end(), paf_cmp);
	std::vector<std::unique_ptr<PAFobject>>::iterator it = paf_objects.begin();
	std::vector<std::unique_ptr<PAFobject>>::iterator next;
	while (it != --paf_objects.end()) {
		next = std::next(it);
		while ((*next)->t_end <= (*it)->t_end) {
			next = paf_objects.erase(next);
			if(next == paf_objects.end()) return;
		}
		it = next;
	}
}

int main(int argc, char** argv) {

	std::vector<std::unique_ptr<PAFobject>> paf_objects;
	auto paf_parser = bioparser::createParser<bioparser::PafParser, PAFobject>(argv[1]);
	paf_parser->parse(paf_objects, -1);

	sweepLineAlgorithm(paf_objects);

	return 0;

}