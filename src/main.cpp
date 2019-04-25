#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <sstream>
#include <tuple>

#include "Config.h"

#include "bioparser/bioparser.hpp"
#include "PAFObject.cpp"
#include "repeats_parser.h"

class FASTAEntity {
	
	public:
		std::string name;
		std::string sequence;
		
		FASTAEntity(
			const char *name, uint32_t name_length,
			const char *sequence, uint32_t sequence_length) {

			(this->name).assign(name, name_length);
			(this->sequence).assign(sequence, sequence_length);
		}
};

void sweepLineAlgorithm(std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
	auto paf_cmp = [](const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) { 
 		if (a->t_begin == b->t_begin) {
        	return (a->t_end > b->t_end);
        }
        return (a->t_begin < b->t_begin);
	};
	std::sort(paf_objects.begin(), paf_objects.end(), paf_cmp);
	std::vector<std::unique_ptr<PAFObject>>::iterator it = paf_objects.begin();
	std::vector<std::unique_ptr<PAFObject>>::iterator next;
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

	std::vector<std::unique_ptr<PAFObject>> paf_objects;
	auto paf_parser = bioparser::createParser<bioparser::PafParser, PAFObject>(argv[1]);
	paf_parser->parse(paf_objects, -1);

	std::vector<std::unique_ptr<FASTAEntity>> ref_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(argv[2]);
	fasta_parser->parse(ref_objects, -1);
	FASTAEntity reference = ref_objects[0];

	sweepLineAlgorithm(paf_objects);

	std::vector<std::tuple<std::string, int, int>> repeats;
	repeats_parser::parse(repeats, argv[3]);

	repeats_parser::remove_covered(repeats, paf_objects);
	std::cout << repeats.size() << std::endl;
	std::tuple<std::string, int, int> r = repeats[0];
	std::cout << std::get<1>(r) << " " << std::get<2>(r) << std::endl;

	return 0;

}