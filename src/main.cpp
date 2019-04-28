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

class Vertex {

public:
	PAFObject read;
	std::vector<Vertex*> vertices;
	Vertex* parent = NULL;

	Vertex(PAFObject read) : read(read) {
		this->read = read;
	}
};

std::vector<Vertex*> create_graph(std::vector<Vertex> &vertices, std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
	for (auto const& paf: paf_objects) {
		vertices.emplace_back(Vertex(*paf));
	}
	std::vector<Vertex>::iterator it = vertices.begin();
	std::vector<Vertex>::iterator next;
	std::vector<Vertex*> heads;
	heads.emplace_back(&(*it));
	while(it != --vertices.end()) {
		next = std::next(it);
		while ((*next).read.t_begin <= (*it).read.t_end) {
			//std::cout << (*it).read.t_begin << std::endl;
			it->vertices.emplace_back(&(*next));
			next++;
			if(next == vertices.end()) break;
		}
		if (next == std::next(it)) heads.emplace_back(&(*next));
		it++;
	}
	FILE* file = fopen ("network.txt","w");
	for (auto const& vertex : vertices) {
		fprintf (file, "%d %d", vertex.read.t_begin, vertex.read.t_end);
		for (auto const& next: vertex.vertices) {
			fprintf (file, " %d", next->read.t_begin);
		}
		fprintf (file, "\n");
	}
	fclose (file);
	return heads;
}

Vertex* DepthFirstSearch(std::vector<Vertex*> heads) {
	Vertex* max;
	Vertex* head;
	Vertex* longest;
	int length=0;
	int begin;
	for (int i = 0; i < heads.size(); i++) {
		head = heads[i];
		begin = head->read.t_begin;
		while(!head->vertices.empty()) {
			max = head->vertices[0];
			for (auto const& vertex : head->vertices) {
				if (vertex->read.t_begin > max->read.t_begin) max = vertex;
			}
			max->parent = head;
			head = max;
		}
		if (length < head->read.t_begin - begin) {
			length = head->read.t_begin - begin;
			longest = head;
		}
	}
	return longest;
}

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

void statistics(Vertex* end, FASTAEntity &reference) {
	int count = 1;
	FILE* file = fopen ("used_reads.txt","w");
	int last_index = end->read.t_end;
	while (end->parent != NULL) {
		count++;
		fprintf (file, "%s %d %d\n", end->read.query_name.c_str(), end->read.t_begin, end->read.t_end);
		end = end->parent;
	}
	fclose (file);
	printf("Genome coverage: %f%\n", (last_index - end->read.t_begin) / (float) reference.sequence.length() * 100);
	printf("Number of used reads: %d\n", count);
}

int main(int argc, char** argv) {

	std::vector<std::unique_ptr<PAFObject>> paf_objects;
	auto paf_parser = bioparser::createParser<bioparser::PafParser, PAFObject>(argv[1]);
	paf_parser->parse(paf_objects, -1);
	sweepLineAlgorithm(paf_objects);

	std::vector<std::unique_ptr<FASTAEntity>> ref_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(argv[2]);
	fasta_parser->parse(ref_objects, -1);
	FASTAEntity reference = *ref_objects[0];

	std::vector<std::tuple<std::string, int, int>> repeats;
	repeats_parser::parse(repeats, argv[3]);
	repeats_parser::remove_covered(repeats, paf_objects);

	if (!repeats_parser::check_repeats(repeats, reference.sequence)) {
		// two same repeats that are not covered by any sequences
		printf("Genome can't be assembled\n");
	}

	std::vector<Vertex> vertices;
	std::vector<Vertex*> heads = create_graph(vertices, paf_objects);

	printf("Number of components: %lu\n", heads.size());

	Vertex* end = DepthFirstSearch(heads);

	statistics(end, reference);

	return 0;

}


