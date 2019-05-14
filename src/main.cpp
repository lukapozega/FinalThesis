#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <sstream>
#include <tuple>
#include <unordered_map>

#include "Config.h"

#include "bioparser/bioparser.hpp"
#include "PAFObject.cpp"
#include "repeats_parser.h"

struct option options[] = {
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

std::unordered_map<char, char> complement_map = {
			{'A', 'T'},
			{'T', 'A'},
			{'G', 'C'},
			{'C', 'G'},
	};

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

bool paf_unique(const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) {
	return a->query_name == b->query_name;
}

void clear_contained_reads(std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
	std::vector<std::unique_ptr<PAFObject>>::iterator it = paf_objects.begin();
	it = std::unique (paf_objects.begin(), paf_objects.end(), paf_unique);
	paf_objects.resize(std::distance(paf_objects.begin(), it));
	paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p){return p->q_end - p->q_begin < 0.9 * p->q_length;}), paf_objects.end());
	auto paf_cmp = [](const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) {
		if (a->target_name == b->target_name) { 
	 		if (a->t_begin == b->t_begin) {
	        	return (a->t_end > b->t_end);
	        }
	        return (a->t_begin < b->t_begin);
	    }
	    return (a->target_name > b->target_name);
	};
	std::sort(paf_objects.begin(), paf_objects.end(), paf_cmp);
	it = paf_objects.begin();
	int s, e;
	std::string n;
	while (it != --paf_objects.end()) {
		s = (*it)->t_begin;
		e = (*it)->t_end;
		n = (*it)->target_name;
		paf_objects.erase(std::remove_if(it+1, paf_objects.end(), [&s, &e, &n](std::unique_ptr<PAFObject> &p){return p->t_begin > s && p->t_end < e && p->target_name == n;}), paf_objects.end());
		printf("%lu\n", paf_objects.size());
		it++;
	}
}

std::string complement(std::string const &sequence){
	std::string comp = "";
	for (char const & c : sequence) {
		comp = complement_map.at(c) + comp; 
	}
	return comp;
}

void statistics(Vertex* end, FASTAEntity &reference, std::string const &filePath, std::string const &reference_name) {
	std::vector<std::unique_ptr<FASTAEntity>> reads;
	auto reads_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(filePath);
	reads_parser->parse(reads, -1);
	int count = 1;
	std::string part;
	int startIndex = end->read.q_begin;
	int endIndex = end->read.q_end;
	FILE* file = fopen ("used_reads.txt","w");
	FILE* genome = fopen ("genome.fasta","w");
	int last_index = end->read.t_end;
	int overlap = 0;
	std::string seq;
	fprintf (genome, ">%s\n", reference_name.c_str());
	while (end->parent != NULL) {
		count++;
		for (auto const& r : reads) {
			if (r->name == end->read.query_name) {
				if (end->read.orientation == '+') {
					seq = complement(r->sequence);
				} else {
					seq = r->sequence;
				}
				part.assign(seq, startIndex + overlap, endIndex - (startIndex + overlap));
			}
		}
		fprintf (genome, "%s", part.c_str());
		endIndex = end->read.q_begin;
		fprintf (file, "%s %d %d\n", end->read.query_name.c_str(), end->read.t_begin, end->read.t_end);
		overlap = end->read.t_begin;
		end = end->parent;
		overlap = end->read.t_end - overlap;
		startIndex = end->read.q_begin;
	}
	for (auto const& r : reads) {
		if (r->name == end->read.query_name) part.assign(r->sequence, startIndex + overlap, endIndex - (startIndex + overlap));
	}
	fprintf (genome, "%s", part.c_str());
	fprintf (file, "%s %d %d\n", end->read.query_name.c_str(), end->read.t_begin, end->read.t_end);
	fclose (file);
	fclose (genome);
	printf("Genome coverage: %f%%\n", (last_index - end->read.t_begin) / (float) reference.sequence.length() * 100);
	printf("Number of used reads: %d\n", count);
}

void help() {
	printf("Program accepts three arguments and prints assemblying statistics.\n");
	printf("Usage: assembly <alignment file file> <reference file> <repeats file> <reads file>\n");
	printf("Arguments should be in PAF, fasta and rpt file format\n");
}

void version() {
	printf("Version %d.%d\n", assembly_VERSION_MAJOR, assembly_VERSION_MINOR);
}

int main(int argc, char** argv) {

	char optchr;
	int option_index = 0;
	while((optchr = getopt_long(argc, argv, "hv", options, &option_index)) != -1) {
		switch (optchr) {
			case 'h': {
				help();
				return 1;
			}
			case 'v': {
				version();
				return 1;
			}
			default: {
				fprintf(stderr, "Unknown option -%c\n", optchr);
				return 1;
			}
		}
	}

	if(argc - optind < 4) {
		fprintf(stderr, "Program requires three arguments.\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	std::vector<std::unique_ptr<PAFObject>> paf_objects;
	auto paf_parser = bioparser::createParser<bioparser::PafParser, PAFObject>(argv[optind]);
	paf_parser->parse(paf_objects, -1);
	clear_contained_reads(paf_objects);

	std::vector<std::unique_ptr<FASTAEntity>> ref_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(argv[optind + 1]);
	fasta_parser->parse(ref_objects, -1);
	FASTAEntity reference = *ref_objects[0];

	std::vector<std::tuple<std::string, int, int>> repeats;
	if (!repeats_parser::parse(repeats, argv[3])) {
		fprintf(stderr, "Error reading file %s\n", argv[optind + 2]);
		return 1;
	}

	repeats_parser::remove_covered(repeats, paf_objects);

	auto repeats_result = repeats_parser::check_repeats(repeats, reference.sequence);

	if (std::get<0>(repeats_result) != -1) {
		printf("Genome can't be assembled\n");
		printf("Non covered areas are matching:\n");
		printf("%d-%d\n", std::get<0>(repeats_result),std::get<1>(repeats_result));
		printf("%d-%d\n", std::get<2>(repeats_result),std::get<3>(repeats_result));
		return 0;
	}

	std::vector<Vertex> vertices;
	std::vector<Vertex*> heads = create_graph(vertices, paf_objects);

	printf("Number of components: %lu\n", heads.size());

	Vertex* end = DepthFirstSearch(heads);

	statistics(end, reference, argv[optind + 3], reference.name);

	return 0;

}


